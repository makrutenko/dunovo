#!/usr/bin/env python3
import argparse
import logging
import os
import sys
import unittest
# Add the utils directory to sys.path so we can import errstats.
script_path = os.path.realpath(__file__)
root_dir = os.path.dirname(os.path.dirname(script_path))
sys.path.append(os.path.join(root_dir, 'utils'))
import errstats

DESCRIPTION = """"""


def make_argparser():
  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  volume = parser.add_mutually_exclusive_group()
  volume.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  volume.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  volume.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')

  #TODO: Allow executing specific tests from the command line.
  #      Probably starts with something like this:
  #      for attr_name in dir(FillInGapQualsTest):
  #        if attr_name.startswith('test_'):
  #          result = FillInGapQualsTest(attr_name)
  #      Or maybe implement a TestRunner?

  unittest.main()


def make_tests(cls, data=(), suite=None):
  """Add tests to a `unittest.TestCase` subclass from a list of input/output data.
  The only requirement is that the `TestCase` subclass implements a `make_test()` classmethod which
  accepts a dict from `data` as `**kwargs`, creates a test function, and returns it.
  `suite` is an optional `unittest.TestSuite` which these tests should be added to.
    Give None to skip this."""
  for datum in data:
    test_function = cls.make_test(**datum)
    setattr(cls, 'test_'+datum['name'], test_function)
  if suite:
    case_suite = unittest.TestLoader().loadTestsFromTestCase(cls)
    suite.addTest(case_suite)


########## utils/errstats.py ##########

def slots_obj_to_dict(o):
  d = {}
  for attr in o.__slots__:
    d[attr] = getattr(o, attr)
  return d

errstatsTests = unittest.TestSuite()

class FillInGapQualsTest(unittest.TestCase):
  @classmethod
  def make_test(cls, input=None, output=None, exception=None, **kwargs):
    def test(self):
      if exception is None:
        result = errstats.fill_in_gap_quals(input)
        self.assertEqual(output, result)
      else:
        self.assertRaises(exception, errstats.fill_in_gap_quals, input)
    return test

make_tests(FillInGapQualsTest, suite=errstatsTests, data=(
    {'name':'nogap',  'input':'HFMHPPHH', 'output':'HFMHPPHH'},
    {'name':'empty',  'input':'',         'output':''},
    {'name':'allsp',  'input':' ', 'exception':ValueError},
    {'name':'simple', 'input':'EIVJX LSFLJ', 'output':'EIVJXOLSFLJ'},
    {'name':'2sp',    'input':'MI  VVIXZ',   'output':'MIPPVVIXZ'},
    {'name':'2gap',   'input':'MI   VVI XZ', 'output':'MIPPPVVISXZ'},
    {'name':'start',  'input':' 04ij0a9ejfa', 'output':'B04ij0a9ejfa'},
    {'name':'end',    'input':'04ij0a9ejfa ', 'output':'04ij0a9ejfae'}
  )
)

class GetAlignmentErrorsTest(unittest.TestCase):
  @classmethod
  def make_test(cls, consensus=None, seq_align=None, qual_align=None, qual_thres=None, errors=None,
                count_indels=None, **kwargs):
    def test(self):
      result = errstats.get_alignment_errors(consensus, seq_align, qual_align, qual_thres,
                                             count_indels=count_indels)
      #TODO: This is a temporary kludge for the drift that happened between errstats.py and these
      #      tests. It now returns Error objects instead of dicts. Ideally I should just make Error
      #      objects comparable with `==`. Also, there's a new `read_coord` key which wasn't in the
      #      test data. I should double-check what those values should be.
      result_dicts = []
      for err_obj in result:
        err_dict = slots_obj_to_dict(err_obj)
        del err_dict['read_coord']
        result_dicts.append(err_dict)
      self.assertEqual(errors, result_dicts)
    return test

#TODO: The 'aln_coord' key used to be 'coord', but I changed it to align with the new return values
#      from errstats.py. I'm pretty sure that's right, but not 100%.
make_tests(GetAlignmentErrorsTest, suite=errstatsTests, data=(
    {'consensus':  'GATTACA', 'qual_thres':0, 'count_indels':True,
     'seq_align': ('GATTACA',
                   'GATTACA'), 'name':'NoErrors',
     'qual_align':('IIIIIII',
                   'IIIIIII'), 'errors':[]},
    {'consensus':  'GATTACA', 'qual_thres':0, 'count_indels':True,
     'seq_align': ('GATTACA',
                   'GATCACA'), 'name':'Snv',
     'qual_align':('IIIIIII',
                   'IIIIIII'), 'errors':[{'alt':'C', 'aln_coord':4, 'seq':1, 'type':'SNV', 'passes':True}]},
    {'consensus':  'GATTACA', 'qual_thres':0, 'count_indels':True,
     'seq_align': ('GATTACA',
                   'GAT-ACA'), 'name':'Del',
     'qual_align':('IIIIIII',
                   'III III'), 'errors':[{'alt':1, 'aln_coord':3, 'seq':1, 'type':'del', 'passes':True}]},
    {'consensus':  'GAT-ACA', 'qual_thres':0, 'count_indels':True,
     'seq_align': ('GAT-ACA',
                   'GATTACA'), 'name':'Ins',
     'qual_align':('IIIIIII',
                   'III III'), 'errors':[{'alt':'T', 'aln_coord':3, 'seq':1, 'type':'ins', 'passes':True}]},
    {'consensus':  'GATTACA', 'qual_thres':0, 'count_indels':False,
     'seq_align': ('GATTAGA',
                   'GAT-ACA'), 'name':'NoCountIndels',
     'qual_align':('IIIIIII',
                   'III III'), 'errors':[{'alt':'G', 'type':'SNV', 'aln_coord':6, 'seq':0, 'passes':True}]},
    {'consensus':  'GATYACA', 'qual_thres':0, 'count_indels':True,
     'seq_align': ('GATTACA',
                   'GATGACA'), 'name':'AmbigSnv',
     'qual_align':('IIIIIII',
                   'IIIIIII'), 'errors':[]},
    {'consensus':  'GATYACA', 'qual_thres':0, 'count_indels':True,
     'seq_align': ('GATTACA',
                   'GAT-ACA'), 'name':'AmbigDel',
     'qual_align':('IIIIIII',
                   'III III'), 'errors':[]},
    {'consensus':  'GATTWCA', 'qual_thres':0, 'count_indels':True,
     'seq_align': ('GATTACA',
                   'GAT-ACA'), 'name':'AmbigAfterDel',
     'qual_align':('IIIIIII',
                   'III III'), 'errors':[{'alt':1, 'type':'del', 'seq':1, 'aln_coord':3, 'passes':True}]},
    {'consensus':  'GATTACA', 'qual_thres':30, 'count_indels':True,
     'seq_align': ('GATTCCA',
                   'GATGACA'), 'name':'QualSnv',
     'qual_align':('IIIIIII',
                   'III+III'), 'errors':[{'alt':'C', 'type':'SNV', 'seq':0, 'aln_coord':5, 'passes':True}]},
    {'consensus':  'GATTACA', 'qual_thres':30, 'count_indels':True,
     'seq_align': ('GATTCCA',
                   'GAT-ACA'), 'name':'QualDel',
     'qual_align':('IIIIIII',
                   'II+ +II'), 'errors':[{'alt':'C', 'type':'SNV', 'seq':0, 'aln_coord':5, 'passes':True}]},
    {'consensus':  'GATTACA', 'qual_thres':30, 'count_indels':True,
     'seq_align': ('GATTACA',
                   'GAT-GCA'), 'name':'QualDelSnv',
     'qual_align':('IIIIIII',
                   'II+ I+I'), 'errors':[{'alt':'G', 'type':'SNV', 'seq':1, 'aln_coord':5, 'passes':True}]},
    {'consensus':  'GATTACA', 'qual_thres':0, 'count_indels':True,
     'seq_align': ('GATTACA',
                   'GAT---A'), 'name':'LongDel',
     'qual_align':('IIIIIII',
                   'III   I'), 'errors':[{'alt':3, 'type':'del', 'seq':1, 'aln_coord':3, 'passes':True}]},
    {'consensus':  'GAT---A', 'qual_thres':0, 'count_indels':True,
     'seq_align': ('GATTACA',
                   'GAT---A'), 'name':'LongIns',
     'qual_align':('IIIIIII',
                   'III   I'), 'errors':[{'alt':'TAC', 'type':'ins', 'seq':0, 'aln_coord':3, 'passes':True}]},
    {'consensus':  'GATTACA', 'qual_thres':30, 'count_indels':True,
     'seq_align': ('GATTACA',
                   'GAT-CCA'), 'name':'SnvAfterDel',
     'qual_align':('IIIIIII',
                   'III III'), 'errors':[{'alt':1, 'type':'del', 'seq':1, 'aln_coord':3, 'passes':True},
                                         {'alt':'C', 'type':'SNV', 'seq':1, 'aln_coord':5, 'passes':True}]},
    {'consensus':  'GATTACA', 'qual_thres':30, 'count_indels':True,
     'seq_align': ('GATTACA',
                   'GAT-CCA'), 'name':'QualAfterDel',
     'qual_align':('IIIIIII',
                   'III +II'), 'errors':[{'alt':1, 'type':'del', 'seq':1, 'aln_coord':3, 'passes':True}]},
    {'consensus':  'GATTAC-', 'qual_thres':30, 'count_indels':True,
     'seq_align': ('GATTAC-',
                   'GATTACA'), 'name':'EndIns',
     'qual_align':('IIIIIII',
                   'IIIIII '), 'errors':[{'alt':'A', 'type':'ins', 'seq':1, 'aln_coord':6, 'passes':False}]},
  )
)


def fail(message):
  logging.critical(message)
  if __name__ == '__main__':
    sys.exit(1)
  else:
    raise Exception('Unrecoverable error')


if __name__ == '__main__':
  sys.exit(main(sys.argv))
