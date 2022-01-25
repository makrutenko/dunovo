#!/usr/bin/env python3
import os
import sys
import argparse
import collections
script_path = os.path.realpath(__file__)
root_dir = os.path.dirname(os.path.dirname(script_path))
sys.path.append(root_dir)
from bfx import getreads

STATS = collections.OrderedDict(
  (
    ('pairs',       'Total pairs of reads in input'),
    ('barcodes',    'Unique barcodes (double-standed families)'),
    ('bar_orders',  'Unique barcode/order combinations (single-stranded families)'),
    ('avg_pairs',   'Average # of read pairs per unique barcode'),
    ('stranded_singletons', 'Single-strand singletons (only 1 read with that barcode/order combination)'),
    ('duplex_singletons',   'Double-strand singletons (only 1 read with that barcode)'),
    ('duplexes',    'Families with both strands present'),
    ('passed_sscs', 'Single-strand families with > {min_reads} reads'),
    ('passed_duplexes',     'Families with both strands with > {min_reads} reads'),
  )
)
USAGE = """$ %(prog)s [options] reads_1.fq reads_2.fq
       $ %(prog)s [options] -f families.tsv"""
DESCRIPTION = """Print statistics on the raw duplex sequencing reads."""
EPILOG = """Warning: This tracks all barcodes in a dict, so it can take a lot of memory. A guideline
is about 200 bytes per (12bp) tag. For example, it took about 800MB for a 10GB, 32 million read
dataset with an average of 4 pairs per barcode."""


def make_argparser():
  parser = argparse.ArgumentParser(usage=USAGE, description=DESCRIPTION, epilog=EPILOG)
  parser.add_argument('infile1', metavar='reads_1.fq', nargs='?',
    help='The first mates in the read pairs.')
  parser.add_argument('infile2', metavar='reads_2.fq', nargs='?',
    help='The second mates in the read pairs.')
  parser.add_argument('-f', '--families', metavar='families.tsv',
    help='The output of make-families.awk.')
  parser.add_argument('-t', '--tag-length', dest='tag_len', type=int, default=12)
  parser.add_argument('-c', '--constant-length', dest='const_len', type=int, default=5)
  parser.add_argument('-m', '--min-reads', type=int, default=3,
    help='The minimum number of reads required in each single-stranded family. Default: '
         '%(default)s')
  parser.add_argument('-v', '--validate', dest='check_ids', action='store_true', default=False,
    help='Check the id\'s of the reads to make sure the correct reads are mated into pairs (the '
         'id\'s of mates must be identical). Default: %(default)s.')
  parser.add_argument('-I', '--no-check-ids', dest='check_ids', action='store_false',
    help='Don\'t check read ids.')
  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  if args.families:
    barcodes = {}
    with open(args.families) as infile:
      barcodes = read_families(infile)
  else:
    if not (args.infile1 and args.infile2):
      fail('Error: You must provide either a --families file or two fastq files.')
    with open(args.infile1) as infileh1:
      with open(args.infile2) as infileh2:
        barcodes = read_fastqs(infileh1, infileh2, tag_len=args.tag_len, check_ids=args.check_ids)

  stats = get_stats(barcodes, stats_key=STATS, min_reads=args.min_reads)
  print_stats(stats, stats_key=STATS)


def read_families(infile):
  barcodes = collections.Counter()
  for line in infile:
    fields = line.split('\t')
    barcode = fields[0]
    order = fields[1]
    barcodes[(barcode, order)] += 1
  return barcodes


def read_fastqs(infileh1, infileh2, tag_len=12, check_ids=False):
  reader1 = getreads.getparser(infileh1, filetype='fastq').parser()
  reader2 = getreads.getparser(infileh2, filetype='fastq').parser()
  barcodes = collections.Counter()
  while True:
    try:
      read1 = next(reader1)
      read2 = next(reader2)
    except StopIteration:
      break
    if check_ids and not read_ids_match(read1.id, read2.id):
      raise getreads.FormatError('Read pair mismatch: "{}" and "{}"'.format(read1.id, read2.id))
    alpha = read1.seq[:tag_len]
    beta  = read2.seq[:tag_len]
    if alpha < beta:
      order = 'ab'
      barcode = alpha + beta
    else:
      order = 'ba'
      barcode = beta + alpha
    barcodes[(barcode, order)] += 1
  return barcodes


def read_ids_match(name1, name2):
  id1 = name1.split()[0]
  id2 = name2.split()[0]
  if id1.endswith('/1'):
    id1 = id1[:-2]
  if id2.endswith('/2'):
    id2 = id2[:-2]
  return id1 == id2


def get_stats(barcodes, stats_key=STATS, min_reads=3):
  stats = {'min_reads':min_reads}
  for stat_name in stats_key:
    stats[stat_name] = 0
  barcode_counts = collections.Counter()
  for (barcode, order), count in barcodes.items():
    stats['pairs'] += count
    barcode_counts[barcode] += count
    if count == 1:
      stats['stranded_singletons'] += 1
    if count >= min_reads:
      stats['passed_sscs'] += 1
    if order == 'ab':
      other_order = 'ba'
    else:
      other_order = 'ab'
    if (barcode, other_order) in barcodes:
      other_count = barcodes[(barcode, other_order)]
      # The loop encounters both strands of each duplex, so to only count each duplex once,
      # only increment on one of the strands (ab).
      if order == 'ab':
        stats['duplexes'] += 1
        if count >= min_reads and other_count >= min_reads:
          stats['passed_duplexes'] += 1
  stats['barcodes'] = len(barcode_counts)
  stats['bar_orders'] = len(barcodes)
  stats['avg_pairs'] = stats['pairs']/stats['barcodes']
  stats['duplex_singletons'] = 0
  for count in barcode_counts.values():
    if count == 1:
      stats['duplex_singletons'] += 1
  return stats


def print_stats(stats, stats_key=STATS):
  for stat_name, stat_description in stats_key.items():
    print(stats[stat_name], stat_description.format(**stats), sep='\t')


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
