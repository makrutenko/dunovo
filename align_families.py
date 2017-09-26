#!/usr/bin/env python
from __future__ import division
import os
import sys
import time
import logging
import tempfile
import argparse
import subprocess
import collections
import distutils.spawn
import seqtools
import parallel
import shims
# There can be problems with the submodules, but none are essential.
# Try to load these modules, but if there's a problem, load a harmless dummy and continue.
simplewrap = shims.get_module_or_shim('utillib.simplewrap')
version = shims.get_module_or_shim('utillib.version')
phone = shims.get_module_or_shim('ET.phone')

#TODO: Warn if it looks like the two input FASTQ files are the same (i.e. the _1 file was given
#      twice). Can tell by whether the alpha and beta (first and last 12bp) portions of the barcodes
#      are always identical. This would be a good thing to warn about, since it's an easy mistake
#      to make, but it's not obvious that it happened. The pipeline won't fail, but will just
#      produce pretty weird results.

REQUIRED_COMMANDS = ['mafft']
USAGE = '$ %(prog)s [options] families.tsv > families.msa.tsv'
DESCRIPTION = """Read in sorted FASTQ data and do multiple sequence alignments of each family."""

def make_argparser():

  wrapper = simplewrap.Wrapper()
  wrap = wrapper.wrap
  parser = argparse.ArgumentParser(usage=USAGE, description=wrap(DESCRIPTION),
                                   formatter_class=argparse.RawTextHelpFormatter)

  wrapper.width = wrapper.width - 24
  parser.add_argument('infile', metavar='read-families.tsv', nargs='?', default=sys.stdin,
                      type=argparse.FileType('r'),
    help=wrap('The input reads, sorted into families. One line per read pair, 8 tab-delimited '
              'columns:\n'
              '1. canonical barcode\n'
              '2. barcode order ("ab" for alpha+beta, "ba" for beta-alpha)\n'
              '3. read 1 name\n'
              '4. read 1 sequence\n'
              '5. read 1 quality scores\n'
              '6. read 2 name\n'
              '7. read 2 sequence\n'
              '8. read 2 quality scores'))
  parser.add_argument('-a', '--aligner', choices=('mafft', 'kalign'), default='mafft',
    help=wrap('The multiple sequence aligner to use.'))
  parser.add_argument('-p', '--processes', type=int, default=1,
    help=wrap('Number of worker subprocesses to use. Must be at least 1. Default: %(default)s.'))
  parser.add_argument('--phone-home', action='store_true',
    help=wrap('Report helpful usage data to the developer, to better understand the use cases and '
              'performance of the tool. The only data which will be recorded is the name and '
              'version of the tool, the size of the input data, the time taken to process it, and '
              'the IP address of the machine running it. No parameters or filenames are sent. All '
              'the reporting and recording code is available at https://github.com/NickSto/ET.'))
  parser.add_argument('--galaxy', dest='platform', action='store_const', const='galaxy',
    help=wrap('Tell the script it\'s running on Galaxy. Currently this only affects data reported '
              'when phoning home.'))
  parser.add_argument('--test', action='store_true',
    help=wrap('If reporting usage data, mark this as a test run.'))
  parser.add_argument('--version', action='version', version=str(version.get_version()),
    help=wrap('Print the version number and exit.'))
  parser.add_argument('-L', '--log-file', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  parser.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
                      default=logging.WARNING)
  parser.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  parser.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)

  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log_file, level=args.volume, format='%(message)s')
  tone_down_logger()

  start_time = time.time()
  if args.phone_home:
    run_id = phone.send_start(__file__, version.get_version(), platform=args.platform,
                              test=args.test, fail='warn')

  assert args.processes > 0, '-p must be greater than zero'

  # Check for required commands.
  missing_commands = []
  for command in REQUIRED_COMMANDS:
    if not distutils.spawn.find_executable(command):
      missing_commands.append(command)
  if missing_commands:
    fail('Error: Missing commands: "'+'", "'.join(missing_commands)+'".')

  # Open a pool of worker processes.
  pool = parallel.StreamingPool(args.processes, process_duplex, [args.aligner])

  # Main loop.
  """This processes whole duplexes (pairs of strands) at a time for a future option to align the
  whole duplex at a time.
  duplex data structure:
  duplex = {
    'ab': [
      {'name1': 'read_name1a',
       'seq1':  'GATT-ACA',
       'qual1': 'sc!0 /J*',
       'name2': 'read_name1b',
       'seq2':  'ACTGACTA',
       'qual2': '34I&SDF)'
      },
      {'name1': 'read_name2a',
       ...
      }
    ]
  }
  e.g.:
  seq = duplex[order][pair_num]['seq1']
  """
  stats = {'duplexes':0, 'time':0, 'pairs':0, 'runs':0, 'failures':0, 'aligned_pairs':0}
  duplex = collections.OrderedDict()
  family = []
  barcode = None
  order = None
  for line in args.infile:
    fields = line.rstrip('\r\n').split('\t')
    if len(fields) != 8:
      continue
    (this_barcode, this_order, name1, seq1, qual1, name2, seq2, qual2) = fields
    # If the barcode or order has changed, we're in a new family.
    # Process the reads we've previously gathered as one family and start a new family.
    if this_barcode != barcode or this_order != order:
      duplex[order] = family
      # If the barcode is different, we're at the end of the whole duplex. Process the it and start
      # a new one. If the barcode is the same, we're in the same duplex, but we've switched strands.
      if this_barcode != barcode:
        # logging.debug('processing {}: {} orders ({})'.format(barcode, len(duplex),
        #               '/'.join([str(len(duplex[o])) for o in duplex])))
        results = pool.compute(duplex, barcode)
        process_results(results, stats)
        duplex = collections.OrderedDict()
      barcode = this_barcode
      order = this_order
      family = []
    pair = {'name1': name1, 'seq1':seq1, 'qual1':qual1, 'name2':name2, 'seq2':seq2, 'qual2':qual2}
    family.append(pair)
    stats['pairs'] += 1
  # Process the last family.
  duplex[order] = family
  # logging.debug('processing {}: {} orders ({}) [last]'.format(barcode, len(duplex),
  #               '/'.join([str(len(duplex[o])) for o in duplex])))
  results = pool.compute(duplex, barcode)
  process_results(results, stats)

  # Process all remaining families in the queue.
  logging.info('flushing..')
  results = pool.flush()
  process_results(results, stats)
  logging.info(pool.states)
  pool.stop()

  if args.infile is not sys.stdin:
    args.infile.close()

  end_time = time.time()
  run_time = int(end_time - start_time)

  # Final stats on the run.
  stats['duplexes'] = pool.jobs_submitted
  logging.error('Processed {pairs} read pairs in {duplexes} duplexes, with {failures} alignment '
                'failures.'.format(**stats))
  if stats['aligned_pairs'] > 0 and stats['runs'] > 0:
    per_pair = stats['time'] / stats['aligned_pairs']
    per_run = stats['time'] / stats['runs']
    logging.error('{:0.3f}s per pair, {:0.3f}s per run.'.format(per_pair, per_run))
  logging.error('in {}s total time.'.format(run_time))

  if args.phone_home:
    stats['align_time'] = stats['time']
    del stats['time']
    phone.send_end(__file__, version.get_version(), run_id, run_time, stats, platform=args.platform,
                   test=args.test, fail='warn')


def process_duplex(duplex, barcode, aligner='mafft'):
  output = ''
  run_stats = {'time':0, 'runs':0, 'aligned_pairs':0, 'failures':0}
  orders = duplex.keys()
  if len(duplex) == 0 or None in duplex:
    return '', {}
  elif len(duplex) == 1:
    # If there's only one strand in the duplex, just process the first mate, then the second.
    combos = ((1, orders[0]), (2, orders[0]))
  elif len(duplex) == 2:
    # If there's two strands, process in a criss-cross order:
    # strand1/mate1, strand2/mate2, strand1/mate2, strand2/mate1
    combos = ((1, orders[0]), (2, orders[1]), (2, orders[0]), (1, orders[1]))
  else:
    raise AssertionError('Error: More than 2 orders in duplex {}: {}'.format(barcode, orders))
  for mate, order in combos:
    family = duplex[order]
    start = time.time()
    try:
      alignment = align_family(family, mate, aligner=aligner)
    except AssertionError as error:
      logging.critical('AssertionError on family {}, order {}, mate {}:\n{}.'
                       .format(barcode, order, mate, error))
      raise
    except (OSError, subprocess.CalledProcessError) as error:
      logging.warning('{} on family {}, order {}, mate {}:\n{}'
                      .format(type(error).__name__, barcode, order, mate, error))
      alignment = None
    # Compile statistics.
    elapsed = time.time() - start
    pairs = len(family)
    logging.info('{} sec for {} read pairs.'.format(elapsed, pairs))
    if pairs > 1:
      run_stats['time'] += elapsed
      run_stats['runs'] += 1
      run_stats['aligned_pairs'] += pairs
    if alignment is None:
      logging.warning('Error aligning family {}/{} (read {}).'.format(barcode, order, mate))
      run_stats['failures'] += 1
    else:
      output += format_msa(alignment, barcode, order, mate)
  return output, run_stats


def align_family(family, mate, aligner='mafft'):
  """Do a multiple sequence alignment of the reads in a family and their quality scores."""
  mate = str(mate)
  assert mate == '1' or mate == '2'
  if len(family) == 0:
    return None
  elif len(family) == 1:
    # If there's only one read pair, there's no alignment to be done (and MAFFT won't accept it).
    aligned_seqs = [family[0]['seq'+mate]]
  else:
    # Do the multiple sequence alignment.
    aligned_seqs = make_msa(family, mate, aligner=aligner)
  # Transfer the alignment to the quality scores.
  ## Get a list of all quality scores in the family for this mate.
  quals_raw = [pair['qual'+mate] for pair in family]
  qual_alignment = seqtools.transfer_gaps_multi(quals_raw, aligned_seqs, gap_char_out=' ')
  # Package them up in the output data structure.
  alignment = []
  for pair, aligned_seq, aligned_qual in zip(family, aligned_seqs, qual_alignment):
    alignment.append({'name':pair['name'+mate], 'seq':aligned_seq, 'qual':aligned_qual})
  return alignment


def make_msa(family, mate, aligner='mafft'):
  if aligner == 'mafft':
    return make_msa_mafft(family, mate)
  elif aligner == 'kalign':
    return make_msa_kalign(family, mate)


def make_msa_kalign(family, mate):
  logging.info('Aligning with kalign.')
  import kalign
  seqs = [pair['seq'+mate] for pair in family]
  return kalign.align(seqs)


def make_msa_mafft(family, mate):
  """Perform a multiple sequence alignment on a set of sequences and parse the result.
  Uses MAFFT."""
  logging.info('Aligning with mafft.')
  #TODO: Replace with tempfile.mkstemp()?
  with tempfile.NamedTemporaryFile('w', delete=False, prefix='align.msa.') as family_file:
    for pair in family:
      name = pair['name'+mate]
      seq = pair['seq'+mate]
      family_file.write('>'+name+'\n')
      family_file.write(seq+'\n')
  with open(os.devnull, 'w') as devnull:
    try:
      command = ['mafft', '--nuc', '--quiet', family_file.name]
      output = subprocess.check_output(command, stderr=devnull)
    except (OSError, subprocess.CalledProcessError):
      raise
    finally:
      # Make sure we delete the temporary file.
      os.remove(family_file.name)
  return read_fasta(output)


def read_fasta(fasta):
  """Quick and dirty FASTA parser. Return the sequences and their names.
  Returns a list of sequences.
  Warning: Reads the entire contents of the file into memory at once."""
  sequences = []
  sequence = ''
  for line in fasta.splitlines():
    if line.startswith('>'):
      if sequence:
        sequences.append(sequence.upper())
      sequence = ''
      continue
    sequence += line.strip()
  if sequence:
    sequences.append(sequence.upper())
  return sequences


def format_msa(align, barcode, order, mate, outfile=sys.stdout):
  output = ''
  for sequence in align:
    output += '{bar}\t{order}\t{mate}\t{name}\t{seq}\t{qual}\n'.format(bar=barcode, order=order,
                                                                       mate=mate, **sequence)
  return output


def process_results(results, stats):
  """Process the outcome of a duplex run.
  Print the aligned output and sum the stats from the run with the running totals."""
  for result in results:
    output, run_stats = result
    for key, value in run_stats.items():
      stats[key] += value
    if output:
      sys.stdout.write(output)


def tone_down_logger():
  """Change the logging level names from all-caps to capitalized lowercase.
  E.g. "WARNING" -> "Warning" (turn down the volume a bit in your log files)"""
  for level in (logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG):
    level_name = logging.getLevelName(level)
    logging.addLevelName(level, level_name.capitalize())


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)


if __name__ == '__main__':
  sys.exit(main(sys.argv))
