#!/usr/bin/env python3
import os
import sys
import gzip
import time
import logging
import argparse
import resource
import subprocess
import networkx
import parallel_tools
import swalign
import shims
from bfx import samreader
from bfx import getreads
# There can be problems with the submodules, but none are essential.
# Try to load these modules, but if there's a problem, load a harmless dummy and continue.
version = shims.get_module_or_shim('utillib.version')
phone = shims.get_module_or_shim('ET.phone')

VERBOSE = (logging.DEBUG+logging.INFO)//2
USAGE = '$ %(prog)s [options] families.tsv barcodes.fa barcodes.sam > families.corrected.tsv'
DESCRIPTION = """Correct barcodes using an alignment of all barcodes to themselves. Reads the
alignment in SAM format and corrects the barcodes in an input "families" file (the output of
make-barcodes.awk). It will print the "families" file to stdout with barcodes (and orders)
corrected."""


def make_argparser():
  parser = argparse.ArgumentParser(usage=USAGE, description=DESCRIPTION)
  parser.add_argument('families', type=open_as_text_or_gzip,
    help="The sorted output of make-barcodes.awk. The important part is that it's a tab-delimited "
      'file with at least 2 columns: the barcode sequence and order, and it must be sorted in the '
      'same order as the "reads" in the SAM file.')
  parser.add_argument('reads', type=open_as_text_or_gzip,
    help='The fasta/q file given to the aligner. Used to get barcode sequences from read names.')
  parser.add_argument('sam', type=argparse.FileType('r'), nargs='?', default=sys.stdin,
    help='Barcode alignment, in SAM format. Omit to read from stdin. The read names must be '
      'integers, representing the (1-based) order they appear in the families file.')
  parser.add_argument('-P', '--prepend', action='store_true',
    help='Prepend the corrected barcodes and orders to the original columns.')
  parser.add_argument('-d', '--dist', type=int, default=1,
    help='NM edit distance threshold. Default: %(default)s')
  parser.add_argument('-m', '--mapq', type=int, default=20,
    help='MAPQ threshold. Default: %(default)s')
  parser.add_argument('-p', '--pos', type=int, default=2,
    help='POS tolerance. Alignments will be ignored if abs(POS - 1) is greater than this value. '
      'Set to greater than the barcode length for no threshold. Default: %(default)s')
  parser.add_argument('-c', '--choose-by', choices=('count', 'connect'), default='count',
    help='Choose the "correct" barcode in a network of related barcodes by either the count of how '
      'many times the barcode was observed ("freq") or how connected the barcode is to the others '
      'in the network ("connect").')
  parser.add_argument('-N', '--allow-no-nm-if-ns', action='store_true',
    help='Allow alignments with missing NM tags if the barcode has at least one N. Otherwise this '
      'will fail if it encounters an alignment missing an NM tag.')
  parser.add_argument('-I', '--no-check-ids', dest='check_ids', action='store_false', default=True,
    help="Don't check to make sure read pairs have identical ids. By default, if this encounters a "
      "pair of reads in families.tsv with ids that aren't identical (minus an ending /1 or /2), it "
      "will throw an error.")
  parser.add_argument('--limit', type=int,
    help='Limit the number of entries that will be read from each input file, for testing purposes.')
  parser.add_argument('-S', '--structures', action='store_true',
    help='Print a list of the unique isoforms')
  parser.add_argument('--struct-human', action='store_true')
  parser.add_argument('-V', '--visualize', metavar='networks.png', nargs='?', default=0,
    help='Produce a visualization of the unique structures and write the image to this file. '
      'If you omit a filename, it will be displayed in a window.')
  parser.add_argument('-F', '--viz-format', choices=('dot', 'graphviz', 'png'), default='png')
  parser.add_argument('-n', '--no-output', dest='output', action='store_false', default=True)
  parser.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  parser.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  parser.add_argument('-i', '--info', dest='volume', action='store_const', const=logging.INFO)
  parser.add_argument('-v', '--verbose', dest='volume', action='store_const', const=VERBOSE)
  parser.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG,
    help='Print debug messages (very verbose).')
  parser.add_argument('--phone-home', action='store_true',
    help='Report helpful usage data to the developer, to better understand the use cases and '
      'performance of the tool. The only data which will be recorded is the name and version of '
      'the tool, the size of the input data, the time and memory taken to process it, and the IP '
      'address of the machine running it. Also, if the script fails, it will report the name of '
      'the exception thrown and the line of code it occurred in. No parameters or filenames are '
      'sent. All the reporting and recording code is available at https://github.com/NickSto/ET.')
  parser.add_argument('--galaxy', dest='platform', action='store_const', const='galaxy',
    help="Tell the script it's running on Galaxy. Currently this only affects data reported when "
      'phoning home.')
  parser.add_argument('--test', action='store_true',
    help='If reporting usage data, mark this as a test run.')
  parser.add_argument('--version', action='version', version=str(version.get_version()),
    help='Print the version number and exit.')
  return parser


def main(argv):

  # Allow using -v for --version if it's the only argument, and --verbose if there are more.
  if len(argv) == 2 and argv[1] == '-v':
    print(version.get_version())
    return

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')
  tone_down_logger()

  start_time = time.time()
  # If the user requested, report back some data about the start of the run.
  if args.phone_home:
    call = phone.Call(
      __file__, version.get_version(), platform=args.platform, test=args.test, fail='warn'
    )
    call.send_data('start')
    call.send_data('prelim', run_data=gather_prelim_data(args.families, args.reads, args.sam))

  # Execute as much of the script as possible in a try/except to catch any exception that occurs
  # and report it via ET.phone.
  try:
    logging.info('Reading the fasta/q to map read names to barcodes..')
    names_to_barcodes = map_names_to_barcodes(args.reads, args.limit)

    logging.info('Reading the SAM to build the graph of barcode relationships..')
    passing_alignments = filter_alignment(
      args.sam, args.pos, args.mapq, args.dist, args.limit, args.allow_no_nm_if_ns
    )
    graph, reversed_barcodes, num_good_alignments = read_alignments(
      passing_alignments, names_to_barcodes
    )

    logging.info('Reading the families.tsv to get the counts of each family..')
    family_counts, read_pairs = get_family_counts(
      args.families, limit=args.limit, check_ids=args.check_ids
    )

    if args.structures or args.visualize != 0:
      logging.info('Counting the unique barcode networks..')
      structures = count_structures(graph, family_counts)
      if args.structures:
        print_structures(structures, args.struct_human)
      if args.visualize != 0:
        logging.info('Generating a visualization of barcode networks..')
        visualize([s['graph'] for s in structures], args.visualize, args.viz_format)

    logging.info('Building the correction table from the graph..')
    corrections = make_correction_table(graph, family_counts, args.choose_by)

    logging.info('Reading the families.tsv again to print corrected output..')
    with open_as_text_or_gzip(args.families.name) as families:
      print_corrected_output(
        families, corrections, reversed_barcodes, args.prepend, args.limit, args.output
      )

    run_time = int(time.time() - start_time)
    max_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
    logging.info(f'Max memory usage: {max_mem:0.2f}MB')
    logging.info(f'Wall clock time:  {run_time} seconds')

  except (Exception, KeyboardInterrupt) as exception:
    if args.phone_home and call:
      try:
        exception_data = getattr(exception, 'child_context', parallel_tools.get_exception_data())
        logging.critical(parallel_tools.format_traceback(exception_data))
        exception_data = parallel_tools.scrub_tb_paths(exception_data, script_path=__file__)
      except Exception:
        exception_data = {}
        raise
      run_time = int(time.time() - start_time)
      try:
        run_data = {
          'barcodes':len(names_to_barcodes), 'good_alignments':num_good_alignments,
          'read_pairs':read_pairs, 'max_mem':int(max_mem)
        }
      except Exception:
        run_data = {}
      run_data['failed'] = True
      if exception_data:
        run_data['exception'] = exception_data
      call.send_data('end', run_time=run_time, run_data=run_data)
      raise exception
    else:
      raise

  if args.phone_home:
    run_data = {
      'barcodes':len(names_to_barcodes), 'good_alignments':num_good_alignments,
      'read_pairs':read_pairs, 'max_mem':int(max_mem)
    }
    call.send_data('end', run_time=run_time, run_data=run_data)


def gather_prelim_data(families, reads, sam):
  data = {}
  data['families_gzipped'] = isinstance(families, gzip.GzipFile)
  data['families_size'] = os.path.getsize(families.name)
  data['reads_gzipped'] = isinstance(reads, gzip.GzipFile)
  data['reads_size'] = os.path.getsize(reads.name)
  data['sam_stdin'] = sam is sys.stdin
  if data['sam_stdin']:
    data['sam_size'] = None
  else:
    data['sam_size'] = os.path.getsize(sam.name)
  return data


def detect_format(reads_file, max_lines=7):
  """Detect whether a file is a fastq or a fasta, based on its content."""
  fasta_votes = 0
  fastq_votes = 0
  line_num = 0
  for line in reads_file:
    line_num += 1
    if line_num % 4 == 1:
      if line.startswith('@'):
        fastq_votes += 1
      elif line.startswith('>'):
        fasta_votes += 1
    elif line_num % 4 == 3:
      if line.startswith('+'):
        fastq_votes += 1
      elif line.startswith('>'):
        fasta_votes += 1
    if line_num >= max_lines:
      break
  reads_file.seek(0)
  if fasta_votes > fastq_votes:
    return 'fasta'
  elif fastq_votes > fasta_votes:
    return 'fastq'
  else:
    return None


def read_fastaq(reads_file):
  filename = reads_file.name.lower()
  if filename.endswith('.fa') or filename.endswith('.fasta'):
    format = 'fasta'
  elif filename.endswith('.fq') or filename.endswith('.fastq'):
    format = 'fastq'
  else:
    format = detect_format(reads_file)
  return getreads.getparser(reads_file, filetype=format)


def map_names_to_barcodes(reads_file, limit=None):
  """Map barcode names to their sequences."""
  names_to_barcodes = {}
  read_num = 0
  for read in read_fastaq(reads_file):
    read_num += 1
    if limit is not None and read_num > limit:
      break
    try:
      name = int(read.name)
    except ValueError:
      logging.critical(f'Non-int read name {read.name!r}')
      raise
    names_to_barcodes[name] = read.seq
  reads_file.close()
  return names_to_barcodes


def filter_alignment(
    sam_file, pos_thres, mapq_thres, dist_thres, limit=None, allow_no_nm_if_ns=False
  ):
  """Read the SAM file and yield reads that pass the filters.
  Returns (qname, rname, reversed)."""
  for aln_num, aln in enumerate(samreader.read(sam_file), 1):
    if limit is not None and aln_num > limit:
      break
    logging.debug(f'read {aln.rname} -> ref {aln.qname} (read seq {aln.seq}):')
    rname_fields = aln.rname.split(':')
    if len(rname_fields) == 2 and rname_fields[1] == 'rev':
      reversed = True
      rname_str = rname_fields[0]
    else:
      reversed = False
      rname_str = aln.rname
    try:
      qname = int(aln.qname)
      rname = int(rname_str)
    except ValueError:
      if aln.rname == '*':
        logging.debug('\tRead unmapped (reference == "*")')
        continue
      else:
        logging.error(
          f'Non-integer read name(s) in alignment {aln_num}: {aln.qname!r}, {rname_str!r}.'
        )
        raise
    if qname == rname:
      logging.debug('\tRead aligned to itself.')
      continue
    # Apply alignment quality filters.
    if aln.flag is None or aln.pos is None or aln.mapq is None:
      logging.warning(f'\tMissing flag ({aln.flag!r}), pos ({aln.pos!r}), or mapq ({aln.mapq!r})')
      continue
    if aln.unmapped:
      logging.debug('\tRead unmapped')
      continue
    if abs(aln.pos - 1) > pos_thres:
      logging.debug(f'\tAlignment failed pos filter: abs({aln.pos} - 1) > {pos_thres}')
      continue
    if aln.mapq < mapq_thres:
      logging.debug(f'\tAlignment failed mapq filter: {aln.mapq} > {mapq_thres}')
      continue
    nm = aln.tags.get('NM')
    if nm is None:
      if 'N' in aln.seq:
        if allow_no_nm_if_ns:
          logging.info(f"\tAlignment missing NM tag, but likely due to N's in sequence: {aln.seq!r}")
          continue
        else:
          raise RuntimeError(f'Alignment missing NM tag in alignment {aln_num}')
    if nm > dist_thres:
      logging.debug(f'\tAlignment failed NM distance filter: {nm} > {dist_thres}')
      continue
    yield qname, rname, reversed


def read_alignments(alignments, names_to_barcodes):
  """Read the alignments from the SAM file.
  Returns (graph, reversed_barcodes, num_good_alignments):
  graph: A networkx.Graph() containing a node per barcode (the sequence as a str), and an edge
    between every pair of barcodes that align to each other (with a threshold-passing alignment).
  reversed_barcodes: The set() of all barcode sequences that are involved in an alignment where the
    target is reversed (swapped halves, like alpha+beta -> beta+alpha). Both the query and reference
    sequence in each alignment are marked here.
  num_good_alignments: The raw number of alignments processed that passed the filters."""
  graph = networkx.Graph()
  reversed_barcodes = set()
  # Maps correct barcode numbers to sets of original barcodes (includes correct ones).
  num_good_alignments = 0
  for qname, rname, reversed in alignments:
    num_good_alignments += 1
    rseq = names_to_barcodes[rname]
    qseq = names_to_barcodes[qname]
    # Is this an alignment to a reversed barcode?
    if reversed:
      reversed_barcodes.add(rseq)
      reversed_barcodes.add(qseq)
    graph.add_node(rseq)
    graph.add_node(qseq)
    graph.add_edge(rseq, qseq)
  return graph, reversed_barcodes, num_good_alignments


def get_family_counts(families_file, limit=None, check_ids=True):
  """For each family (barcode), count how many read pairs exist for each strand (order)."""
  family_counts = {}
  last_barcode = None
  this_family_counts = None
  read_pairs = 0
  for line in families_file:
    read_pairs += 1
    if limit is not None and read_pairs > limit:
      break
    fields = line.rstrip('\r\n').split('\t')
    if check_ids:
      assert_read_ids_match(fields[2], fields[5])
    barcode = fields[0]
    order = fields[1]
    if barcode != last_barcode:
      if this_family_counts:
        this_family_counts['all'] = this_family_counts['ab'] + this_family_counts['ba']
      family_counts[last_barcode] = this_family_counts
      this_family_counts = {'ab':0, 'ba':0}
      last_barcode = barcode
    this_family_counts[order] += 1
  this_family_counts['all'] = this_family_counts['ab'] + this_family_counts['ba']
  family_counts[last_barcode] = this_family_counts
  families_file.close()
  return family_counts, read_pairs


def assert_read_ids_match(name1, name2):
  id1 = name1.split()[0]
  id2 = name2.split()[0]
  if id1.endswith('/1'):
    id1 = id1[:-2]
  if id2.endswith('/2'):
    id2 = id2[:-2]
  if id1 == id2:
    return True
  elif id1.endswith('/2') and id2.endswith('/1'):
    raise ValueError(
      f'Read names not as expected. Mate 1 ends with /2 and mate 2 ends with /1:\n'
      f'  Mate 1: {name1!r}\n  Mate 2: {name2!r}'
    )
  else:
    raise ValueError(f'Read names {name1!r} and {name2!r} do not match.')


def make_correction_table(meta_graph, family_counts, choose_by='count'):
  """Make a table mapping original barcode sequences to correct barcodes.
  Assumes the most connected node in the graph as the correct barcode."""
  corrections = {}
  for nodes in networkx.connected_components(meta_graph):
    graph = meta_graph.subgraph(nodes)
    if choose_by == 'count':
      def key(bar):
        return family_counts[bar]['all']
    elif choose_by == 'connect':
      degrees = graph.degree()
      def key(bar):
        return degrees[bar]
    barcodes = sorted(graph.nodes(), key=key, reverse=True)
    correct = barcodes[0]
    for barcode in barcodes:
      if barcode != correct:
        logging.debug(f'Correcting {barcode} ->\n           {correct}\n')
        corrections[barcode] = correct
  return corrections


def print_corrected_output(
    families_file, corrections, reversed_barcodes, prepend=False, limit=None, output=True
  ):
  line_num = 0
  barcode_num = 0
  barcode_last = None
  corrected = {'reads':0, 'barcodes':0, 'reversed':0}
  reads = [0, 0]
  corrections_in_this_family = 0
  for line in families_file:
    line_num += 1
    if limit is not None and line_num > limit:
      break
    fields = line.rstrip('\r\n').split('\t')
    raw_barcode = fields[0]
    order = fields[1]
    if raw_barcode != barcode_last:
      # We just started a new family.
      barcode_num += 1
      family_info = f'{barcode_last}\t{reads[0]}\t{reads[1]}'
      if corrections_in_this_family:
        corrected['reads'] += corrections_in_this_family
        corrected['barcodes'] += 1
        family_info += '\tCORRECTED!'
      else:
        family_info += '\tuncorrected'
      logging.log(VERBOSE, family_info)
      reads = [0, 0]
      corrections_in_this_family = 0
      barcode_last = raw_barcode
    if order == 'ab':
      reads[0] += 1
    elif order == 'ba':
      reads[1] += 1
    if raw_barcode in corrections:
      correct_barcode = corrections[raw_barcode]
      corrections_in_this_family += 1
      # Check if the order of the barcode reverses in the correct version.
      # First, we check in reversed_barcodes whether either barcode was involved in a reversed
      # alignment, to save time (is_alignment_reversed() does a full smith-waterman alignment).
      if ((raw_barcode in reversed_barcodes or correct_barcode in reversed_barcodes) and
          is_alignment_reversed(raw_barcode, correct_barcode)):
        # If so, then switch the order field.
        corrected['reversed'] += 1
        if order == 'ab':
          correct_order = 'ba'
        else:
          correct_order = 'ab'
      else:
        correct_order = order
    else:
      correct_barcode = raw_barcode
      correct_order = order
    # Add the corrected barcode and order to the output.
    if prepend:
      fields.insert(0, correct_barcode)
      fields.insert(1, correct_order)
    else:
      fields[0] = correct_barcode
      fields[1] = correct_order
    if output:
      print(*fields, sep='\t')
  if corrections_in_this_family:
    corrected['reads'] += corrections_in_this_family
    corrected['barcodes'] += 1
  logging.info(
    'Corrected {barcodes} barcodes on {reads} read pairs, with {reversed} reversed.'
    .format(**corrected)
  )


def is_alignment_reversed(barcode1, barcode2):
  """Return True if the barcodes are reversed with respect to each other, False otherwise.
  "reversed" in this case meaning the alpha + beta halves are swapped.
  Determine by aligning the two to each other, once in their original forms, and once with the
  second barcode reversed. If the smith-waterman score is higher in the reversed form, return True.
  """
  half = len(barcode2)//2
  barcode2_rev = barcode2[half:] + barcode2[:half]
  fwd_align = swalign.smith_waterman(barcode1, barcode2)
  rev_align = swalign.smith_waterman(barcode1, barcode2_rev)
  if rev_align.score > fwd_align.score:
    return True
  else:
    return False


def count_structures(meta_graph, family_counts):
  """Count the number of unique (isomorphic) subgraphs in the main graph."""
  structures = []
  for nodes in networkx.connected_components(meta_graph):
    graph = meta_graph.subgraph(nodes)
    match = False
    for structure in structures:
      archetype = structure['graph']
      if networkx.is_isomorphic(graph, archetype):
        match = True
        structure['count'] += 1
        structure['central'] += int(is_centralized(graph, family_counts))
        break
    if not match:
      size = len(graph)
      central = is_centralized(graph, family_counts)
      structures.append({'graph':graph, 'size':size, 'count':1, 'central':int(central)})
  return structures


def is_centralized(graph, family_counts):
  """Checks if the graph is centralized in terms of where the reads are located.
  In a centralized graph, the node with the highest degree is the only one which (may) have more
  than one read pair associated with that barcode.
  This returns True if that's the case, False otherwise."""
  if len(graph) == 2:
    # Special-case graphs with 2 nodes, since the other algorithm doesn't work for them.
    # - When both nodes have a degree of 1, sorting by degree doesn't work and can result in the
    #   barcode with more read pairs coming second.
    barcode1, barcode2 = graph.nodes()
    counts1 = family_counts[barcode1]
    counts2 = family_counts[barcode2]
    total1 = counts1['all']
    total2 = counts2['all']
    logging.debug(
      f'{barcode1}: {total1:3d} ({counts1["ab"]}/{counts1["ba"]})\n'
      f'{barcode2}: {total2:3d} ({counts2["ab"]}/{counts2["ba"]})'
    )
    if (total1 >= 1 and total2 == 1) or (total1 == 1 and total2 >= 1):
      return True
    else:
      return False
  else:
    degrees = graph.degree()
    first = True
    for barcode in sorted(graph.nodes(), key=lambda barcode: degrees[barcode], reverse=True):
      if not first:
        counts = family_counts[barcode]
        # How many read pairs are associated with this barcode (how many times did we see this barcode)?
        try:
          if counts['all'] > 1:
            return False
        except TypeError:
          logging.critical(f'barcode: {barcode}, counts: {counts}')
          raise
      first = False
    return True


def print_structures(structures, human=True):
  # Define a cmp function to sort the list of structures in ascending order of size, but then
  # descending order of count.
  def cmp_fxn(structure1, structure2):
    if structure1['size'] == structure2['size']:
      return structure2['count'] - structure1['count']
    else:
      return structure1['size'] - structure2['size']
  width = None
  last_size = None
  for structure in sorted(structures, cmp=cmp_fxn):
    size = structure['size']
    graph = structure['graph']
    if size == last_size:
      i += 1
    else:
      i = 0
    if width is None:
      width = str(len(str(structure['count'])))
    letters = num_to_letters(i)
    degrees = sorted(graph.degree().values(), reverse=True)
    if human:
      degrees_str = ' '.join(map(str, degrees))
    else:
      degrees_str = ','.join(map(str, degrees))
    if human:
      format_str = '{:2d}{:<3s} {count:<'+width+'d} {central:<'+width+'d} {}'
      print(format_str.format(size, letters+':', degrees_str, **structure))
    else:
      print(size, letters, structure['count'], structure['central'], degrees_str, sep='\t')
    last_size = size


def num_to_letters(i):
  """Translate numbers to letters, e.g. 1 -> A, 10 -> J, 100 -> CV"""
  letters = ''
  while i > 0:
    n = (i-1) % 26
    i = i // 26
    if n == 25:
      i -= 1
    letters = chr(65+n) + letters
  return letters


def visualize(graphs, viz_path, args_viz_format):
  import matplotlib
  from networkx.drawing.nx_agraph import graphviz_layout
  meta_graph = networkx.Graph()
  for graph in graphs:
    add_graph(meta_graph, graph)
  pos = graphviz_layout(meta_graph)
  networkx.draw(meta_graph, pos)
  if viz_path:
    ext = os.path.splitext(viz_path)[1]
    if ext == '.dot':
      viz_format = 'graphviz'
    elif ext == '.png':
      viz_format = 'png'
  else:
    viz_format = args_viz_format
  if viz_format == 'graphviz':
    from networkx.drawing.nx_pydot import write_dot
    assert viz_path is not None, 'Must provide a filename to --visualize if using --viz-format "graphviz".'
    base_path = os.path.splitext(viz_path)
    write_dot(meta_graph, base_path+'.dot')
    run_command('dot', '-T', 'png', '-o', base_path+'.png', base_path+'.dot')
    logging.info('Wrote image of graph to '+base_path+'.dot')
  elif viz_format == 'png':
    if viz_path is None:
      matplotlib.pyplot.show()
    else:
      matplotlib.pyplot.savefig(viz_path)


def add_graph(graph, subgraph):
  # I'm sure there's a function in the library for this, but just cause I need it quick..
  for node in subgraph.nodes():
    graph.add_node(node)
  for edge in subgraph.edges():
    graph.add_edge(*edge)
  return graph


def open_as_text_or_gzip(path):
  """Return an open file-like object reading the path as a text file or a gzip file, depending on
  which it looks like."""
  if detect_gzip(path):
    return gzip.open(path, 'r')
  else:
    return open(path, 'r')


def detect_gzip(path):
  """Return True if the file looks like a gzip file: ends with .gz or contains non-ASCII bytes."""
  ext = os.path.splitext(path)[1]
  if ext == '.gz':
    return True
  elif ext in ('.txt', '.tsv', '.csv'):
    return False
  with open(path) as fh:
    is_not_ascii = detect_non_ascii(fh.read(100))
  if is_not_ascii:
    return True


def detect_non_ascii(bytes, max_test=100):
  """Return True if any of the first "max_test" bytes are non-ASCII (the high bit set to 1).
  Return False otherwise."""
  for i, char in enumerate(bytes):
    # Is the high bit a 1?
    if ord(char) & 128:
      return True
    if i >= max_test:
      return False
  return False


def run_command(*command):
  try:
    exit_status = subprocess.call(command)
  except subprocess.CalledProcessError as cpe:
    exit_status = cpe.returncode
  except OSError:
    exit_status = None
  return exit_status


def tone_down_logger():
  """Change the logging level names from all-caps to capitalized lowercase.
  E.g. "WARNING" -> "Warning" (turn down the volume a bit in your log files)"""
  for level in (logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG):
    level_name = logging.getLevelName(level)
    logging.addLevelName(level, level_name.capitalize())


if __name__ == '__main__':
  sys.exit(main(sys.argv))
