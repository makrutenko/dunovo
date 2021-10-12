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

# A pair of `ReadFamily`s with the same order and barcode.
class StrandFamily(typing.NamedTuple):
  order: str
  mate1: ReadFamily
  mate2: ReadFamily

# A pair of `StrandFamily`s with the same barcode.
class BarFamily(typing.NamedTuple):
  bar: str
  ab: StrandFamily
  ba: StrandFamily


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
    read1s.append(Read(name1, seq1, quals1))
    read2s.append(Read(name2, seq2, quals2))
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
