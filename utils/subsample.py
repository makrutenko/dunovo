#!/usr/bin/env python3
import argparse
import logging
import pathlib
import random
import sys
script_dir = pathlib.Path(__file__).resolve().parent
root_dir = script_dir.parent
sys.path.append(str(root_dir))
import dunovo_parsers
from bfx import getreads

USAGE = "%(prog)s [options]"
DESCRIPTION = """"""


def make_argparser():
  parser = argparse.ArgumentParser(add_help=False, description=DESCRIPTION)
  options = parser.add_argument_group('Options')
  parser.add_argument('families', metavar='families.tsv', type=argparse.FileType('r'),
    help='The input reads, sorted into families (the output of make-families.sh).')
  parser.add_argument('fastq1', metavar='reads_1.fq', type=argparse.FileType('r'),
    help='The raw input reads (mate 1).')
  parser.add_argument('fastq2', metavar='reads_2.fq', type=argparse.FileType('r'),
    help='The raw input reads (mate 2).')
  parser.add_argument('fastq1_out', metavar='selected_reads_1.fq', type=argparse.FileType('w'),
    help='Output filename for the subsampled reads (mate 1). '
      'Warning: any existing file will be overwritten.')
  parser.add_argument('fastq2_out', metavar='selected_reads_2.fq', type=argparse.FileType('w'),
    help='Output filename for the subsampled reads (mate 2). '
      'Warning: any existing file will be overwritten.')
  parser.add_argument('-f', '--fraction', type=float, default=1,
    help='Fraction of families to output. Default: %(default)s')
  parser.add_argument('-s', '--seed', type=int,
    help='Random number generator seed. Default: %(default)s')
  parser.add_argument('-P', '--prepended', action='store_true',
    help='The families.tsv file is the result of correct.py --prepend.')
  options.add_argument('-h', '--help', action='help',
    help='Print this argument help text and exit.')
  logs = parser.add_argument_group('Logging')
  logs.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  volume = logs.add_mutually_exclusive_group()
  volume.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  volume.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  volume.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')

  if args.seed is None:
    seed = random.randint(1, 2**32)
    logging.error(f'Seed: {seed}')
  else:
    seed = args.seed
  random.seed(seed)

  logging.info('Info: Reading input families file..')
  families = dunovo_parsers.parse_make_families(args.families, args.prepended)
  chosen_families = choose_elements(families, args.fraction)
  mate1_names, mate2_names = get_all_read_names(chosen_families, args.fraction)
  logging.info('Info: Filtering raw FASTQ into subsampled FASTQ (mate 1)..')
  find_and_write_chosen_reads(mate1_names, args.fastq1, args.fastq1_out)
  logging.info('Info: Filtering raw FASTQ into subsampled FASTQ (mate 2)..')
  find_and_write_chosen_reads(mate2_names, args.fastq2, args.fastq2_out)


def choose_elements(elements, fraction):
  for element in elements:
    if random.random() <= fraction:
      yield element


def get_all_read_names(families, fraction):
  all_mate1_names = set()
  all_mate2_names = set()
  for family in families:
    mate1_names, mate2_names = get_read_names(family)
    all_mate1_names.update(mate1_names)
    all_mate2_names.update(mate2_names)
  return all_mate1_names, all_mate2_names


def get_read_names(family):
  mate1_names = []
  mate2_names = []
  for strand_family in family.ab, family.ba:
    mate1s = strand_family.mate1.reads
    mate1_names.extend([read.name for read in mate1s])
    mate2s = strand_family.mate2.reads
    mate2_names.extend([read.name for read in mate2s])
  return mate1_names, mate2_names


def find_and_write_chosen_reads(chosen_names, input_fastq, output_fastq):
  input_reads = getreads.getparser(input_fastq, filetype='fastq')
  chosen_reads = find_chosen_reads(input_reads, chosen_names)
  write_reads(chosen_reads, output_fastq)


def find_chosen_reads(input_reads, read_names):
  for read in input_reads:
    if read.name in read_names:
      yield read


def write_reads(reads, outfile):
  for read in reads:
    print(read.to_fastq(), file=outfile)


def fail(message):
  logging.critical('Error: '+str(message))
  if __name__ == '__main__':
    sys.exit(1)
  else:
    raise Exception(message)


if __name__ == '__main__':
  try:
    sys.exit(main(sys.argv))
  except BrokenPipeError:
    pass
