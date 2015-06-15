#!/usr/bin/env python
from __future__ import division
import os
import sys
import time
import logging
import tempfile
import argparse
import subprocess
import distutils.spawn

REQUIRED_COMMANDS = ['mafft', 'em_cons']
OPT_DEFAULTS = {'min_reads':3, 'processes':1}
USAGE = "%(prog)s [options]"
DESCRIPTION = """Build single-strand consensus sequences from read families. Pipe sorted reads into
stdin. Prints single-strand consensus sequences to stdout. The sequence names are BARCODE.MATE, e.g.
"CTCAGATAACATACCTTATATGCA.1"."""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('infile', metavar='read-families.tsv', nargs='?',
    help='The input reads, sorted into families.')
  parser.add_argument('-r', '--min-reads', type=int,
    help='The minimum number of reads required to form a family. Families with fewer reads will '
         'be skipped. Default: %(default)s.')
  parser.add_argument('-s', '--stats-file',
    help='Print statistics on the run to this file. Use "-" to print to stderr.')
  parser.add_argument('-p', '--processes', type=int,
    help='Number of processes to use. If > 1, launches this many worker subprocesses. '
         'Default: %(default)s.')
  parser.add_argument('-S', '--slurm', action='store_true',
    help='If -p > 1, prepend sub-commands with "srun -C new".')

  args = parser.parse_args(argv[1:])

  assert args.processes > 0

  # Check for required commands.
  missing_commands = []
  if args.slurm:
    REQUIRED_COMMANDS.append('srun')
  for command in REQUIRED_COMMANDS:
    if not distutils.spawn.find_executable(command):
      missing_commands.append(command)
  if missing_commands:
    fail('Error: Missing commands: "'+'", "'.join(missing_commands)+'".')

  if args.infile:
    infile = open(args.infile)
  else:
    infile = sys.stdin

  if args.stats_file:
    if args.stats_file == '-':
      logging.basicConfig(stream=sys.stderr, level=logging.INFO, format='%(message)s')
    else:
      logging.basicConfig(filename=args.stats_file, filemode='w', level=logging.INFO,
                          format='%(message)s')
  else:
    logging.disable(logging.CRITICAL)

  # Open all the worker processes, if we're using more than one.
  if args.processes > 1:
    workers = open_workers(args.processes, slurm=args.slurm, stats_file=args.stats_file)

  total_time = 0
  total_pairs = 0
  total_runs = 0
  all_pairs = 0
  all_families = 0
  family = []
  family_barcode = None
  for line in infile:
    fields = line.rstrip('\r\n').split('\t')
    if len(fields) != 7:
      continue
    (barcode, name1, seq1, qual1, name2, seq2, qual2) = fields
    # If the barcode has changed, we're in a new family.
    # Process the reads we've previously gathered as one family and start a new family.
    if barcode != family_barcode:
      if len(family) >= args.min_reads:
        all_families += 1
        if args.processes == 1:
          (elapsed, pairs) = process_family(family, family_barcode)
          if pairs > 1:
            total_time += elapsed
            total_pairs += pairs
            total_runs += 1
        else:
          i = all_families % len(workers)
          worker = workers[i]
          delegate(worker, family, family_barcode)
      family_barcode = barcode
      family = []
    family.append((name1, seq1, qual1, name2, seq2, qual2))
    all_pairs += 1
  # Process the last family.
  if len(family) >= args.min_reads:
    all_families += 1
    if args.processes == 1:
      (elapsed, pairs) = process_family(family, family_barcode)
      if pairs > 1:
        total_time += elapsed
        total_pairs += pairs
        total_runs += 1
    else:
      i = all_families % len(workers)
      worker = workers[i]
      delegate(worker, family, family_barcode)

  if args.processes > 1:
    close_workers(workers)
    compile_results(workers)
    # delete_tempfiles(workers)

  if infile is not sys.stdin:
    infile.close()

  if not args.stats_file:
    return

  # Final stats on the run.
  logging.info('Processed {} read pairs and {} multi-pair families.'.format(all_pairs, total_runs))
  per_pair = total_time / total_pairs
  per_run = total_time / total_runs
  logging.info('{:0.3f}s per pair, {:0.3f}s per run.'.format(per_pair, per_run))


def open_workers(num_workers, slurm=False, stats_file=None):
  """Open the required number of worker processes."""
  script_path = os.path.realpath(sys.argv[0])
  workers = []
  for i in range(num_workers):
    if slurm:
      command = ['srun', '-C', 'new', 'python', script_path]
    else:
      command = ['python', script_path]
    stats_subfile = None
    if stats_file:
      if stats_file == '-':
        stats_subfile = '-'
      else:
        stats_subfile = "{}.{}.log".format(stats_file, i)
      command.extend(['-s', stats_subfile])
    outfile = open("tmp-sscs-out.{}.fa".format(i), 'w')
    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=outfile)
    worker = {'proc':process, 'outfile':outfile, 'stats':stats_subfile}
    workers.append(worker)
  return workers


def delegate(worker, family, barcode):
  """Send a family to a worker process."""
  for pair in family:
    line = barcode+'\t'+'\t'.join(pair)+'\n'
    worker['proc'].stdin.write(line)


def close_workers(workers):
  for worker in workers:
    worker['outfile'].close()
    worker['proc'].stdin.close()


def compile_results(workers):
  for worker in workers:
    with open(worker['outfile'].name, 'r') as outfile:
      for line in outfile:
        sys.stdout.write(line)


def delete_tempfiles(workers):
  for worker in workers:
    os.remove(worker['outfile'].name)
    if worker['stats']:
      os.remove(worker['stats'])


def process_family(family, barcode):
  start = time.time()
  pairs = len(family)
  if pairs == 1:
    (name1, seq1, qual1, name2, seq2, qual2) = family[0]
    print '>'+barcode+'.1'
    print seq1
    print '>'+barcode+'.2'
    print seq2
  else:
    align_path = make_msa(family, 1)
    consensus = get_consensus(align_path)
    if consensus is not None:
      print '>'+barcode+'.1'
      print consensus
    align_path = make_msa(family, 2)
    consensus = get_consensus(align_path)
    if consensus is not None:
      print '>'+barcode+'.2'
      print consensus
  end = time.time()
  elapsed = end - start
  logging.info('{} sec for {} read pairs.'.format(elapsed, pairs))
  return (elapsed, pairs)


def make_msa(family, mate):
  """Perform a multiple sequence alignment on a set of sequences.
  Uses MAFFT."""
  #TODO: Replace with tempfile.mkstemp()?
  with tempfile.NamedTemporaryFile('w', delete=False, prefix='sscs.') as family_file:
    for pair in family:
      if mate == 1:
        name = pair[0]
        seq = pair[1]
      else:
        name = pair[3]
        seq = pair[4]
      family_file.write('>'+name+'\n')
      family_file.write(seq+'\n')
  with tempfile.NamedTemporaryFile('w', delete=False, prefix='sscs.') as align_file:
    with open(os.devnull, 'w') as devnull:
      command = ['mafft', '--nuc', '--quiet', family_file.name]
      subprocess.call(command, stdout=align_file, stderr=devnull)
  os.remove(family_file.name)
  return align_file.name


def get_consensus(align_path):
  """Make a consensus from a multiple sequence alignment file and return the
  consensus sequence as a string.
  Uses the EMBOSS em_cons command."""
  # Note on em_cons output:
  # It may always be lowercase, but maybe not. It can contain "N", and possibly "-".
  with tempfile.NamedTemporaryFile('w', delete=False, prefix='sscs.') as cons_file:
    cons_path = cons_file.name
  with open(os.devnull, 'w') as devnull:
    command = ['em_cons', '-sequence', align_path, '-outseq', cons_path]
    subprocess.call(command, stderr=devnull)
  os.remove(align_path)
  if os.path.getsize(cons_path) == 0:
    os.remove(cons_path)
    return None
  else:
    consensus = read_fasta(cons_path)
    os.remove(cons_path)
    return consensus


def read_fasta(fasta_path):
  """Read a FASTA file, return the sequence.
  Uses a very narrow definition of FASTA: That returned by the "em_cons" command."""
  seq_lines = []
  at_header = True
  with open(fasta_path) as fasta_file:
    for line in fasta_file:
      if at_header and line.startswith('>'):
        at_header = False
        continue
      seq_lines.append(line.strip())
  return "".join(seq_lines)


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
