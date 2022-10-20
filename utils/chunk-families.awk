# Usage: awk -v nchunks=4 -v outbase=path/to/outfiles -f chunk-families.awk families.tsv families.tsv

BEGIN {
  FS="\t"
  OFS="\t"
  chunk_num = 1
  if (! outbase) {
    outbase = "family-chunks"
  }
  outfile = outbase "." chunk_num ".tsv"
  printf("") > outfile
}

NR > FNR {
  if (! total_lines) {
    total_lines = NR - FNR
    chunk_size = int(total_lines / nchunks)
    chunk_line = FNR - 1
    end_chunk = 0
  }
  chunk_line ++
  if (chunk_line == chunk_size) {
    end_chunk = 1
    final_family = $1
  }
  if (end_chunk && $1 != final_family) {
    chunk_num ++
    outfile = outbase "." chunk_num ".tsv"
    printf("") > outfile
    chunk_line = 1
    end_chunk = 0
  }
  print >> outfile
}