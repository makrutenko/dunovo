#!/usr/bin/env python3
import os
import sys
import time
import logging
import argparse
import resource
import collections
import parallel_tools
import consensus
import swalign
import shims
# There can be problems with the submodules, but none are essential.
# Try to load these modules, but if there's a problem, load a harmless dummy and continue.
simplewrap = shims.get_module_or_shim('utillib.simplewrap')
version = shims.get_module_or_shim('utillib.version')
phone = shims.get_module_or_shim('ET.phone')

# The ascii values that represent a 0 PHRED score.
QUAL_OFFSETS = {'sanger':33, 'solexa':64}
USAGE = """$ %(prog)s [options] families.msa.tsv -1 duplexes_1.fa -2 duplexes_2.fa
       $ cat families.msa.tsv | %(prog)s [options] -1 duplexes_1.fa -2 duplexes_2.fa"""
DESCRIPTION = """Build consensus sequences from read aligned families. Prints duplex consensus \
sequences in FASTA to stdout. The sequence ids are BARCODE.MATE, e.g. "CTCAGATAACATACCTTATATGCA.1", \
where "BARCODE" is the input barcode, and "MATE" is "1" or "2" as an arbitrary designation of the \
two reads in the pair. The id is followed by the count of the number of reads in the two families \
(one from each strand) that make up the duplex, in the format READS1/READS2. If the duplex is \
actually a single-strand consensus because the matching strand is missing, only one number is \
listed.
Rules for consensus building: Single-strand consensus sequences are made by counting how many of \
each base are at a given position. Bases with a PHRED quality score below the --qual threshold are \
not counted. If a majority of the reads (that pass the --qual threshold at that position) have one \
base at that position, then that base is used as the consensus base. If no base has a majority, then \
an N is used. Duplex consensus sequences are made by aligning pairs of single-strand consensuses, \
and comparing bases at each position. If they agree, that base is used in the consensus. Otherwise, \
the IUPAC ambiguity code for both bases is used (N + anything and gap + non-gap result in an N)."""


def make_argparser():

  wrapper = simplewrap.Wrapper()
  wrap = wrapper.wrap
  parser = argparse.ArgumentParser(usage=USAGE, description=wrap(DESCRIPTION), add_help=False,
                                   formatter_class=argparse.RawTextHelpFormatter)

  wrapper.width = wrapper.width - 24
  io = parser.add_argument_group('Inputs and outputs')
  io.add_argument('infile', metavar='families.msa.tsv', nargs='?', default=sys.stdin,
    type=argparse.FileType('r'),
    help=wrap('The output of align_families.py. 6 columns:\n'
              '1. (canonical) barcode\n'
              '2. order ("ab" or "ba")\n'
              '3. mate ("1" or "2")\n'
              '4. read name\n'
              '5. aligned sequence\n'
              '6. aligned quality scores.'))
  io.add_argument('-1', '--dcs1', metavar='duplex_1.fa', type=argparse.FileType('w'),
    help=wrap('The file to output the first mates of the duplex consensus sequences into. '
              'Warning: This will be overwritten if it exists!'))
  io.add_argument('-2', '--dcs2', metavar='duplex_2.fa', type=argparse.FileType('w'),
    help=wrap('Same, but for mate 2.'))
  io.add_argument('--sscs1', metavar='sscs_1.fa', type=argparse.FileType('w'),
    help=wrap('Save the single-strand consensus sequences (mate 1) in this file (FASTA format). '
              'Warning: This will be overwritten if it exists!'))
  io.add_argument('--sscs2', metavar='sscs_2.fa', type=argparse.FileType('w'),
    help=wrap('Save the single-strand consensus sequences (mate 2) in this file (FASTA format). '
              'Warning: This will be overwritten if it exists!'))
  io.add_argument('-F', '--qual-format', choices=QUAL_OFFSETS.keys(), default='sanger',
    help=wrap('FASTQ quality score format. Sanger scores are assumed to begin at '
              f"{QUAL_OFFSETS['sanger']} ({chr(QUAL_OFFSETS['sanger'])!r}). Default: %(default)s."))
  io.add_argument('--fastq-out', metavar='PHRED_SCORE', type=int,
    help=wrap('Output in FASTQ instead of FASTA. You must specify the quality score to give to all '
              'bases. There is no meaningful quality score we can automatically give, so you will '
              'have to specify an artificial one. A good choice is 40, the maximum score normally '
              'output by sequencers.'))
  params = parser.add_argument_group('Algorithm parameters')
  params.add_argument('-r', '--min-reads', type=int, default=3,
    help=wrap('The minimum number of reads (from each strand) required to form a single-strand '
              'consensus. Strands with fewer reads will be skipped. Default: %(default)s.'))
  params.add_argument('-q', '--qual', metavar='PHRED_SCORE', type=int, default=20,
    help=wrap('Base quality threshold. Bases below this quality will not be counted. '
              'Default: %(default)s.'))
  params.add_argument('-c', '--cons-thres', metavar='THRES', type=float, default=0.5,
    help=wrap('The fractional threshold to use when making consensus sequences. The consensus base '
              'must be present in more than this fraction of the reads, or N will be used as the '
              'consensus base instead. Default: %(default)s'))
  params.add_argument('-C', '--min-cons-reads', type=int, default=0,
    help=wrap('The absolute threshold to use when making consensus sequences. The consensus base '
              'must be present in more than this number of reads, or N will be used as the '
              'consensus base instead. Default: %(default)s'))
  phoning = parser.add_argument_group('Feedback')
  phoning.add_argument('--phone-home', action='store_true',
    help=wrap('Report helpful usage data to the developer, to better understand the use cases and '
              'performance of the tool. The only data which will be recorded is the name and '
              'version of the tool, the size of the input data, the time and memory taken to '
              'process it, and the IP address of the machine running it. Also, if the script '
              'fails, it will report the name of the exception thrown and the line of code it '
              'occurred in. No filenames are sent, and the only parameters reported are the number '
              'of --processes and the --queue-size. All the reporting and recording code is '
              'available at https://github.com/NickSto/ET.'))
  phoning.add_argument('--galaxy', dest='platform', action='store_const', const='galaxy',
    help=wrap('Tell the script it\'s running on Galaxy. Currently this only affects data reported '
              'when phoning home.'))
  phoning.add_argument('--test', action='store_true',
    help=wrap('If reporting usage data, mark this as a test run.'))
  log = parser.add_argument_group('Logging')
  log.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help=wrap('Print log messages to this file instead of to stderr. Warning: Will overwrite the '
              'file.'))
  log.add_argument('-Q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
                   default=logging.WARNING)
  log.add_argument('-V', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  log.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  misc = parser.add_argument_group('Miscellaneous')
  misc.add_argument('-p', '--processes', default=0,
    help=wrap('Number of worker subprocesses to use. If 0, no subprocesses will be started and '
              'everything will be done inside one process. Give "auto" to use as many processes '
              'as there are CPU cores. Default: %(default)s.'))
  misc.add_argument('--queue-size', type=int,
    help=wrap('How long to go accumulating responses from worker subprocesses before dealing '
              f'with all of them. Default: {parallel_tools.QUEUE_SIZE_MULTIPLIER} * the number of '
              'worker --processes.'))
  misc.add_argument('-v', '--version', action='version', version=str(version.get_version()),
    help=wrap('Print the version number and exit.'))
  misc.add_argument('-h', '--help', action='store_true',
    help='Print this text on usage and arguments.')

  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])
  if args.help:
    parser.print_help()
    return 0

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')
  tone_down_logger()

  start_time = time.time()
  # If the user requested, report back some data about the start of the run.
  if args.phone_home:
    call = phone.Call(
      __file__, version.get_version(), platform=args.platform, test=args.test, fail='warn'
    )
    call.send_data('start')
    data = {
      'stdin': args.infile is sys.stdin, 'processes': args.processes, 'queue_size': args.queue_size,
    }
    if data['stdin']:
      data['input_size'] = None
    else:
      data['input_size'] = os.path.getsize(args.infile.name)
    call.send_data('prelim', run_data=data)
  else:
    call = None

  # Execute as much of the script as possible in a try/except to catch any exception that occurs
  # and report it via ET.phone.
  try:
    # Process and validate arguments.
    if args.queue_size is not None and args.queue_size <= 0:
      fail('Error: --queue-size must be greater than zero.')
    qual_start = QUAL_OFFSETS[args.qual_format]
    qual_thres = chr(args.qual + qual_start)
    if args.fastq_out is None:
      # Output FASTA.
      output_qual = None
    else:
      # Output FASTQ.
      if qual_start+args.fastq_out > 126:
        fail(f'Error: --fastq-out PHRED score ({args.fastq_out}) is too large.')
      output_qual = chr(qual_start+args.fastq_out)
    if args.min_cons_reads > args.min_reads:
      fail(
        "Error: --min-reads must be greater than --min-cons-reads (or you'll have a lot of "
        "consensus sequences with only N's!). If you want to exclude families with fewer than X "
        "reads, give --min-reads X instead of --min-cons-reads X."
      )
    if not any((args.dcs1, args.dcs2, args.sscs1, args.sscs2)):
      fail('Error: must specify an output file!')
    # A dict of output filehandles.
    # Indexed so we can do filehandles['dcs'][mate].
    filehandles = {'dcs': (args.dcs1, args.dcs2), 'sscs': (args.sscs1, args.sscs2)}

    # Open a pool of worker processes.
    stats = {'time':0, 'reads':0, 'runs':0, 'duplexes':0}
    static_kwargs = {
      'min_reads': args.min_reads, 'cons_thres': args.cons_thres,
      'min_cons_reads': args.min_cons_reads, 'qual_thres': qual_thres, 'output_qual': output_qual,
    }
    pool = parallel_tools.SyncAsyncPool(
      process_duplex, processes=args.processes, static_kwargs=static_kwargs,
      queue_size=args.queue_size, callback=process_result, callback_args=[filehandles, stats],
    )
    try:
      process_families(args.infile, pool, stats)
    finally:
      # If the root process encounters an exception and doesn't tell the workers to stop, it will
      # hang forever.
      pool.close()
      pool.join()
      # Close all open filehandles.
      if args.infile is not sys.stdin:
        args.infile.close()
      for fh_group in filehandles.values():
        for fh in fh_group:
          if fh:
            fh.close()

    # Final stats on the run.
    run_time = int(time.time() - start_time)
    max_mem = get_max_mem()
    logging.info(
      f"Processed {stats['total_reads']} reads and {stats['runs']} duplexes in {run_time} seconds."
    )
    if stats['reads'] > 0 and stats['runs'] > 0:
      per_read = stats['time'] / stats['reads']
      per_run = stats['time'] / stats['runs']
      logging.info(f'{per_read:0.3f}s per read, {per_run:0.3f}s per run.')
    logging.info(f'in {run_time}s total time and {max_mem:0.2f}MB RAM.')

  except (Exception, KeyboardInterrupt) as exception:
    if args.phone_home and call:
      try:
        exception_data = getattr(exception, 'child_context', parallel_tools.get_exception_data())
        logging.critical(parallel_tools.format_traceback(exception_data))
        exception_data = parallel_tools.scrub_tb_paths(exception_data, script_path=__file__)
      except Exception:
        exception_data = {}
      run_time = int(time.time() - start_time)
      try:
        run_data = get_run_data(stats, pool)
      except (Exception, UnboundLocalError):
        run_data = {}
      try:
        run_data['mem'] = get_max_mem()
      except Exception:
        pass
      run_data['failed'] = True
      if exception_data:
        run_data['exception'] = exception_data
      call.send_data('end', run_time=run_time, run_data=run_data)
      raise exception
    else:
      raise

  if args.phone_home and call:
    run_data = get_run_data(stats, pool, max_mem)
    call.send_data('end', run_time=run_time, run_data=run_data)


def process_families(infile, pool, stats):
  total_reads = 0
  duplex = collections.OrderedDict()
  family = []
  barcode = None
  order = None
  # Note: mate is a 0-indexed integer ("mate 1" from the input file is mate 0 here).
  mate = None
  for line in infile:
    # Allow comments (e.g. for test input files).
    if line.startswith('#'):
      continue
    fields = line.rstrip('\r\n').split('\t')
    if len(fields) != 6:
      continue
    this_barcode, this_order, this_mate_str, name, seq, qual = fields
    this_mate = int(this_mate_str)-1
    # If the barcode, order, and mate are the same, we're just continuing the add reads to the
    # current family. Otherwise, store the current family, start a new one, and process the
    # duplex if we're at the end of one.
    new_barcode = this_barcode != barcode
    new_order = this_order != order
    new_mate = this_mate != mate
    if new_barcode or new_order or new_mate:
      if order is not None and mate is not None:
        duplex[(order, mate)] = family
      # If the barcode changed, process the last duplex and start a new one.
      if new_barcode and barcode is not None:
        assert len(duplex) <= 4, duplex.keys()
        pool.compute(duplex, barcode)
        stats['duplexes'] += 1
        duplex = collections.OrderedDict()
      barcode = this_barcode
      order = this_order
      mate = this_mate
      family = []
    read = {'name': name, 'seq':seq, 'qual':qual}
    family.append(read)
    total_reads += 1
  # Process the last family.
  if order is not None and mate is not None:
    duplex[(order, mate)] = family
  assert len(duplex) <= 4, duplex.keys()
  pool.compute(duplex, barcode)
  stats['duplexes'] += 1
  stats['total_reads'] = total_reads
  # Retrieve the remaining results.
  logging.info('Flushing remaining results from worker processes..')
  pool.flush()


def get_max_mem():
  """Get the maximum memory usage (RSS) of this process and all its children, in MB."""
  maxrss_total  = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
  maxrss_total += resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
  return maxrss_total/1024


def get_run_data(stats, pool, max_mem=None):
  run_data = stats.copy()
  run_data['consensus_time'] = run_data['time']
  del run_data['time']
  if max_mem is not None:
    run_data['mem'] = max_mem
  run_data['processes'] = pool.processes
  run_data['queue_size'] = pool.queue_size
  return run_data


def process_duplex(
  duplex, barcode, min_reads=3, cons_thres=0.5, min_cons_reads=0, qual_thres=' ', output_qual=None
):
  """Create duplex consensus sequences for the reads from one barcode."""
  # The code in the main loop used to ensure that "duplex" contains only reads belonging to one final
  # duplex consensus read: ab.1 and ba.2 reads OR ab.2 and ba.1 reads. (Of course, one half might
  # be missing).
  famsizes = ', '.join([str(len(family)) for family in duplex.values()])
  logging.info(f'Starting duplex {barcode}: {famsizes}')
  start = time.time()
  # Construct consensus sequences.
  try:
    sscss = make_sscss(duplex, min_reads, cons_thres, min_cons_reads, qual_thres)
    dcss = make_dcss(sscss)
  except AssertionError:
    logging.exception(f'While processing duplex {barcode}:')
    raise
  # Format output.
  dcs_strs, sscs_strs = format_outputs(dcss, sscss, barcode, output_qual)
  # Calculate run statistics.
  elapsed = time.time() - start
  total_reads = sum([len(family) for family in duplex.values()])
  logging.debug(f'{elapsed} sec for {total_reads} reads.')
  if len(sscss) > 0:
    run_stats = {'time':elapsed, 'runs':1, 'reads':total_reads}
  else:
    run_stats = {'time':0, 'runs':0, 'reads':0}
  return dcs_strs, sscs_strs, run_stats


def make_sscss(duplex, min_reads, cons_thres, min_cons_reads, qual_thres):
  """Create single-strand consensus sequences from families of raw reads."""
  sscss = {}
  for (order, mate), family in duplex.items():
    # logging.info(f'\t{order}.{mate}:')
    # for read in family:
    #   logging.info('\t\t{name}\t{seq}'.format(**read))
    if len(family) < min_reads:
      logging.debug(f'\tnot enough reads ({len(family)} < {min_reads})')
      continue
    sscs = make_sscs(family, order, mate, qual_thres, cons_thres, min_cons_reads)
    sscss[(order, mate)] = sscs
  return sscss


def make_sscs(family, order, mate, qual_thres, cons_thres, min_cons_reads):
  seqs = [read['seq'] for read in family]
  quals = [read['qual'] for read in family]
  consensus_seq = consensus.get_consensus(
    seqs, quals, cons_thres=cons_thres, min_reads=min_cons_reads, qual_thres=qual_thres
  )
  return {'seq':consensus_seq, 'order':order, 'mate':mate, 'nreads':len(family)}


def make_dcss(sscss):
  # ordermates is the mapping between the duplex consensus mate number and the order/mates of the
  # SSCSs it's composed of. It's arbitrary but consistent, to make sure the duplex consensuses have
  # different mate numbers, and they're the same from run to run.
  ordermates = {
    0: (('ab', 0), ('ba', 1)),
    1: (('ab', 1), ('ba', 0)),
  }
  # Get the consensus of each pair of SSCSs.
  # dcss is indexed by (0-based) mate.
  dcss = []
  for duplex_mate in 0, 1:
    # Gather the pair of reads for this duplex consensus.
    sscs_pair = []
    for order, mate in ordermates[duplex_mate]:
      sscs = sscss.get((order, mate))
      if sscs:
        sscs_pair.append(sscs)
    if len(sscs_pair) < 2:
      # If we didn't find two SSCSs for this duplex mate, we can't make a complete pair of duplex
      # consensus sequences.
      break
    align = swalign.smith_waterman(sscs_pair[0]['seq'], sscs_pair[1]['seq'])
    if len(align.target) != len(align.query):
      message = f'{len(align.target)} != {len(align.query)}:\n'
      message += '\n'.join([repr(sscs) for sscs in sscs_pair])
      raise AssertionError(message)
    seq = consensus.build_consensus_duplex_simple(align.target, align.query)
    reads_per_strand = [sscs['nreads'] for sscs in sscs_pair]
    dcss.append({'seq':seq, 'nreads':reads_per_strand})
  assert len(dcss) == 0 or len(dcss) == 2, len(dcss)
  return dcss


def format_outputs(dcss, sscss, barcode, output_qual=None):
  """Format the consensus sequences into FASTA/Q-formatted strings ready for printing.
  sscs_strs is structured so that sscs_strs[order][mate] is the FASTA/Q-formatted output string for
  one SSCS (including ending newline).
  If output_qual is not None, this will format into FASTQ, using the output_qual as the quality
  score character for every base. Note: output_qual must be a string of length 1."""
  # SSCS
  sscs_strs = {}
  for order in 'ab', 'ba':
    sscs_str_pair = []
    for mate in 0, 1:
      sscs = sscss.get((order, mate))
      if sscs and sscs['seq']:
        if output_qual is None:
          sscs_str_pair.append('>{bar}.{order} {nreads}\n{seq}\n'.format(bar=barcode, **sscs))
        else:
          quals = output_qual * len(sscs['seq'])
          sscs_str_pair.append(
            '@{bar}.{order} {nreads}\n{seq}\n+\n{quals}\n'.format(bar=barcode, quals=quals, **sscs)
          )
    if len(sscs_str_pair) == 2:
      sscs_strs[order] = sscs_str_pair
  # DCS
  dcs_strs = []
  if dcss and all([dcs['seq'] for dcs in dcss]):
    assert len(dcss) == 2, (barcode, len(dcss))
    for duplex_mate in 0, 1:
      dcs = dcss[duplex_mate]
      nreads_str = '-'.join([str(nreads) for nreads in dcs['nreads']])
      if output_qual is None:
        dcs_str = f">{barcode} {nreads_str}\n{dcs['seq']}\n"
      else:
        quals = output_qual * len(dcs['seq'])
        dcs_str = f"@{barcode} {nreads_str}\n{dcs['seq']}\n+\n{quals}\n"
      dcs_strs.append(dcs_str)
  return dcs_strs, sscs_strs


def process_result(result, filehandles, stats):
  dcs_strs, sscs_strs, run_stats = result
  # Stats
  for key, value in run_stats.items():
    stats[key] += value
  # SSCSs
  for order, sscs_pair in sscs_strs.items():
    if sscs_pair:
      for mate in 0, 1:
        sscs_fh = filehandles['sscs'][mate]
        if sscs_fh:
          sscs_fh.write(sscs_pair[mate])
  # DCSs
  if dcs_strs:
    for duplex_mate in 0, 1:
      dcs_fh = filehandles['dcs'][duplex_mate]
      if dcs_fh:
        dcs_fh.write(dcs_strs[duplex_mate])


def tone_down_logger():
  """Change the logging level names from all-caps to capitalized lowercase.
  E.g. "WARNING" -> "Warning" (turn down the volume a bit in your log files)"""
  for level in (logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG):
    level_name = logging.getLevelName(level)
    logging.addLevelName(level, level_name.capitalize())


def fail(message):
  logging.critical(message)
  sys.exit(1)


if __name__ == '__main__':
  sys.exit(main(sys.argv))
