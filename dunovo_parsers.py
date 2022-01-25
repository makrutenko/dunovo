import typing
from bfx.getreads import Read

# NamedTuple has been chosen over dataclasses for performance reasons.
# This module is intended to be usable by the main Du Novo scripts, for which performance and memory
# usage is an important constraint. And it seems dataclasses still have higher memory usage than
# NamedTuples:
# https://stackoverflow.com/questions/51671699/data-classes-vs-typing-namedtuple-primary-use-cases/51673969#51673969

# A set of `Read`s with the same mate, order, and barcode.
class ReadFamily(typing.NamedTuple):
  mate: int
  reads: typing.Iterable[Read]
  def get_summary(self):
    return f'mate{self.mate}: {len(self.reads)}'
  def __getitem__(self, index):
    return self.reads[index]
  def __bool__(self):
    return bool(self.reads)

# A pair of `ReadFamily`s with the same order and barcode.
class StrandFamily(typing.NamedTuple):
  order: str
  mate1: ReadFamily
  mate2: ReadFamily
  def get_summary(self):
    summaries = []
    for mate in 'mate1', 'mate2':
      read_family = getattr(self, mate, None)
      if read_family:
        summaries.append(read_family.get_summary())
      else:
        summaries.append(f'{mate}: None')
    return f'{self.order}: ({summaries[0]}, {summaries[1]})'
  def __getitem__(self, keys):
    try:
      len(keys)
    except TypeError:
      mate = keys
      index = None
    else:
      mate, index = keys
    read_family = getattr(self, f'mate{mate}')
    if read_family is None or index is None:
      return read_family
    return read_family[index]
  def __bool__(self):
    return bool(self.mate1 or self.mate2)

# A pair of `StrandFamily`s with the same barcode.
class BarFamily(typing.NamedTuple):
  bar: str
  ab: StrandFamily
  ba: StrandFamily
  def __getitem__(self, keys):
    mate = index = None
    if isinstance(keys, str):
      order = keys
    else:
      order = keys[0]
      if len(keys) >= 2:
        mate = keys[1]
      if len(keys) >= 3:
        index = keys[2]
    strand_family = getattr(self, order)
    if strand_family is None or mate is None:
      return strand_family
    return strand_family[mate,index]
  def get_summary(self):
    summaries = []
    for order in 'ab', 'ba':
      strand_family = getattr(self, order, None)
      if strand_family:
        summaries.append(strand_family.get_summary())
      else:
        summaries.append(f'{order}: None')
    return f'{self.bar}: ({summaries[0]}, {summaries[1]})'
  def __bool__(self):
    return bool(self.ab or self.ba)


class DunovoFormatError(ValueError):
  pass


def parse_make_families(lines, prepended=False):
  strand_families = []
  strand_family_lines = []
  last_barcode = last_order = None
  if prepended:
    expected_columns = 10
  else:
    expected_columns = 8
  for line_num, line in enumerate(lines, 1):
    fields = line.rstrip('\r\n').split('\t')
    if len(fields) != expected_columns:
      raise DunovoFormatError(f'Line {line_num} has an invalid number of columns: {len(fields)}')
    # If it's the output of correct.py with --prepend, there's an extra column.
    # We want the corrected barcode (column 1), not the original one (column 2).
    if prepended:
      fields[2:4] = []
    barcode, order = fields[:2]
    if barcode != last_barcode or order != last_order:
      if last_order is not None:
        strand_families.append(create_strand_family(strand_family_lines))
      strand_family_lines = []
    if barcode != last_barcode:
      if last_barcode is not None:
        yield create_bar_family(strand_families, last_barcode)
      strand_families = []
    strand_family_lines.append(fields)
    last_barcode = barcode
    last_order = order
  if last_order is not None:
    strand_families.append(create_strand_family(strand_family_lines))
  if last_barcode is not None:
    yield create_bar_family(strand_families, last_barcode)


def create_strand_family(strand_family_lines):
  read1s = []
  read2s = []
  last_order = None
  for fields in strand_family_lines:
    barcode, order, name1, seq1, quals1, name2, seq2, quals2 = fields
    if order not in ('ab', 'ba'):
      raise DunovoFormatError(f'Invalid order: {order!r}')
    assert order == last_order or last_order is None, (order, last_order)
    read1s.append(Read(name=name1, seq=seq1, qual=quals1))
    read2s.append(Read(name=name2, seq=seq2, qual=quals2))
    last_order = order
  read_family1 = ReadFamily(1, tuple(read1s))
  read_family2 = ReadFamily(2, tuple(read2s))
  return StrandFamily(order, read_family1, read_family2)


def create_bar_family(strand_families_raw, barcode):
  assert 1 <= len(strand_families_raw) <= 2, len(strand_families_raw)
  # Create a strand_families list with them in the right order.
  strand_families = [None, None]
  for strand_family in strand_families_raw:
    if strand_family.order == 'ab':
      strand_families[0] = strand_family
    elif strand_family.order == 'ba':
      strand_families[1] = strand_family
  # Fill in any missing strand families with empty ones.
  for i, (order, strand_family) in enumerate(zip(('ab', 'ba'), strand_families)):
    if strand_family is None:
      strand_families[i] = StrandFamily(order, ReadFamily(1,()), ReadFamily(2,()))
  return BarFamily(barcode, *strand_families)


def parse_msa(lines):
  bar_family = strand_family = read_family = None
  last_mate = last_order = last_bar = None
  for line_num, line in enumerate(lines,1):
    # Parse the values from the line.
    fields = line.rstrip('\r\n').split('\t')
    if len(fields) != 6:
      raise DunovoFormatError(f'Line {line_num} has an invalid number of columns: {len(fields)}')
    barcode, order, mate_str, name, seq, quals = fields
    try:
      mate = int(mate_str)
    except ValueError:
      raise DunovoFormatError(f'Line {line_num} has an invalid mate column: {mate_str!r}') from None
    print(f'Parsed line {line_num}: {barcode} {order} {mate}')
    # Reset any families if we've started a new one, and yield the BarFamily if we've finished one.
    if barcode != last_bar:
      if bar_family:
        yield bar_family
      read_family = strand_family = bar_family = None
    if order != last_order:
      read_family = strand_family = None
    if mate != last_mate:
      read_family = None
    # Find the right families for this line, if they exist.
    if bar_family is not None:
      strand_family = getattr(bar_family, order)
    if strand_family is not None:
      read_family = getattr(strand_family, f'mate{mate}')
    # Create any families that don't exist yet.
    if read_family is None:
      print('  read_family was None. Creating.')
      read_family = ReadFamily(mate=mate, reads=[])
    if strand_family is None:
      print('  strand_family was None. Creating.')
      attrs = {'order':order, f'mate{mate}':read_family, f'mate{other_mate(mate)}':None}
      strand_family = StrandFamily(**attrs)
    if bar_family is None:
      print('  bar_family was None. Creating.')
      attrs = {'bar':barcode, order:strand_family, other_order(order):None}
      bar_family = BarFamily(**attrs)
    # Add the read to the ReadFamily and properly nest the families.
    read = Read(name=name, seq=seq, qual=quals)
    read_family.reads.append(read)
    attrs = {f'mate{mate}':read_family}
    strand_family = strand_family._replace(**attrs)
    attrs = {order:strand_family}
    bar_family = bar_family._replace(**attrs)
    # Set the last values to the current values.
    last_bar = barcode
    last_order = order
    last_mate = mate
  if bar_family:
    yield bar_family


def other_mate(mate):
  """Return the number of the other mate.
  E.g. `other_mate(1) == 2` and `other_mate(2) == 1`"""
  return ((mate-1) ^ 1) + 1


def other_order(order):
  if order == 'ab':
    return 'ba'
  elif order == 'ba':
    return 'ab'
