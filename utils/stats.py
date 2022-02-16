#!/usr/bin/env python3
import os
import sys
import math
import argparse
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import seqtools
from bfx import swalign

INF = float('inf')
STATS = ('diffs', 'diffs-binned', 'seqlen', 'strand')
OPT_DEFAULTS = {'bins':10, 'probes':'', 'thres':0.75}
USAGE = "%(prog)s [options]"
DESCRIPTION = """"""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('stats',
    help='The type of statistics to compute and print. Give a comma-separated list of stat names, '
         'choosing from "{}".'.format('", "'.join(STATS)))
  parser.add_argument('infile', metavar='read-families.msa.tsv', nargs='?',
    help='The --msa output of sscs.py. Will read from stdin if not provided.')
  parser.add_argument('-b', '--bins', type=int,
    help='The number of bins to segment reads into when doing "diffs-binned".')
  parser.add_argument('-p', '--probes',
    help='Sequence excerpts from the sense strand. Required for "strand" statistic. '
         'Comma-separated.')
  parser.add_argument('-t', '--thres', type=int,
    help='Alignment identity threshold (in fraction, not decimal). Default: %(default)s')

  args = parser.parse_args(argv[1:])

  stats = args.stats.split(',')
  for stat in stats:
    if stat not in STATS:
      fail('Error: invalid statistic "{}". Must choose one of "{}".'.format(stat, '", "'.join(STATS)))
  if 'strand' in stats and not args.probes:
    fail('Error: must provide a probe if requesting "strand" statistic.')

  if args.infile:
    infile = open(args.infile)
  else:
    infile = sys.stdin

  family = []
  consensus = None
  barcode = None
  for line in infile:
    fields = line.rstrip('\r\n').split('\t')
    if len(fields) != 3:
      continue
    (this_barcode, name, seq) = fields
    if fields[1] == 'CONSENSUS':
      if family and consensus:
        process_family(stats, barcode, consensus, family, args)
      barcode = this_barcode
      consensus = seq
      family = []
    else:
      family.append(seq)
  if family and consensus:
    process_family(stats, barcode, consensus, family, args)

  if infile is not sys.stdin:
    infile.close()


#TODO: Maybe print the number of N's in the consensus?
def process_family(stats, barcode, consensus, family, args):
  # Compute stats requiring the whole family at once.
  for stat in stats:
    if stat == 'diffs':
      diffs = seqtools.get_diffs_frac_simple(consensus, family)
    elif stat == 'diffs-binned':
      diffs_binned = seqtools.get_diffs_frac_binned(consensus, family, args.bins)
    elif stat == 'strand':
      probes = args.probes.split(',')
      strand = get_strand(consensus, probes, args.thres)
  # Print the requested stats for each read.
  # Columns: barcode, [stat columns], read sequence.
  for (i, read) in enumerate(family):
    sys.stdout.write(barcode+'\t')
    for stat in stats:
      if stat == 'diffs':
        sys.stdout.write('{}\t'.format(round_sig_figs(diffs[i], 3)))
      elif stat == 'diffs-binned':
        if diffs_binned is None:
          sys.stdout.write('\t' * args.bins)
        else:
          for diff in diffs_binned[i]:
            sys.stdout.write(str(round_sig_figs(diff, 3))+'\t')
      elif stat == 'seqlen':
        sys.stdout.write('{}\t'.format(len(read)))
      elif stat == 'strand':
        sys.stdout.write('{}\t'.format(strand))
    print(read.upper())


def get_strand(seq, probes, thres):
  """Determine which strand the sequence comes from by trying to align probes from the sense strand.
  Returns 'sense', 'anti', or None.
  Algorithm: This tries each probe in both directions.
  If at least one of the alignments has an identity above the threshold, a vote is cast for the
  direction with a higher identity.
  If the votes that were cast are unanimous for one direction, that strand is returned.
  Else, return None."""
  votes = []
  for probe in probes:
    alignment = swalign.smith_waterman(seq, probe)
    sense_id = alignment.matches/len(probe)
    alignment = swalign.smith_waterman(seq, seqtools.get_revcomp(probe))
    anti_id  = alignment.matches/len(probe)
    # print '{}: sense: {}, anti: {}'.format(probe, sense_id, anti_id)
    if sense_id > thres or anti_id > thres:
      if sense_id > anti_id:
        votes.append('sense')
      else:
        votes.append('anti')
  strand = None
  for vote in votes:
    if strand:
      if strand != vote:
        return None
    else:
      strand = vote
  return strand


def round_sig_figs(n, figs):
  if n == 0:
    return n
  elif n < 0:
    n = -n
    sign = -1
  elif n > 0:
    sign = 1
  elif math.isnan(n) or n == INF:
    return n
  magnitude = int(math.floor(math.log10(n)))
  return sign * round(n, figs - 1 - magnitude)


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)


if __name__ == '__main__':
  sys.exit(main(sys.argv))
