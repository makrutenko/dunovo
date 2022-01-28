#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
# get the name of the test directory
dirname=$(dirname $0)

cmd_prefix="$dirname/../"
Usage="Usage: \$ $(basename $0) [options] [test1 [test2]]
Options:
-p: Assume script executables are on your \$PATH.
    E.g. execute make-consensi.py instead of ${cmd_prefix}make-consensi.py.
-q: Quiet mode.
-v: Verbose mode."

function main {

  gave_tests=
  verbose=true
  # Run the requested tests
  for arg in "$@"; do
    # Check for options
    #TODO: option to keep test data at end instead of removing it.
    if [[ ${arg:0:1} == '-' ]]; then
      case "$arg" in
        -h)
          print_usage >&2
          exit 1;;
        -p)
          cmd_prefix=;;
        -q)
          verbose='';;
        -v)
          verbose=true;;
        *)
          echo "Unrecognized option \"$arg\"." >&2;;
      esac
      continue
    fi
    # Execute valid tests (if they're existing functions).
    if [[ $(type -t $arg) == function ]]; then
      gave_tests=true
      if [[ $verbose ]]; then
        $arg
      else
        $arg 2>/dev/null
      fi
    else
      echo "Unrecognized test \"$arg\"." >&2
    fi
  done

  # If no tests were specified in arguments, do all tests.
  if ! [[ $gave_tests ]]; then
    fail "Error: Please specify a valid test to run. Use -h option to list them."
  fi
}

function print_usage {
  echo "$Usage"
  echo "Meta tests:"
  list_meta_tests | indent
  echo "Active tests:"
  list_active_tests | indent
  echo "Inactive tests:"
  list_inactive_tests | indent
  echo "Core tests:"
  list_core_tests | indent
  echo "Unit tests:"
  echo "  unit"
}

function fail {
  echo "$@" >&2
  exit 1
}

function indent {
  awk '{print "  " $0}'
}

function list_functions {
  declare -F | awk 'substr($3,1,1) != "_" {print $3}'
}

function list_base_tests {
  # Print a list of all tests, excluding meta tests.
  list_functions | while read function; do
    if ! printf '%s\n' "$initial_functions_plus_meta" | grep -qxF "$function"; then
      # Filter out regular functions and meta tests.
      printf '%s\n' "$function"
    fi
  done
}

function list_active_tests {
  list_base_tests | while read function; do
    if printf '%s\n' "$all_functions_minus_inactive" | grep -qxF "$function"; then
      printf '%s\n' "$function"
    fi
  done
}

function list_inactive_tests {
  list_base_tests | while read function; do
    if ! printf '%s\n' "$all_functions_minus_inactive" | grep -qxF "$function"; then
      printf '%s\n' "$function"
    fi
  done
}

function list_core_tests {
  list_base_tests | while read function; do
    # Get functions up to and including the core ones.
    if printf '%s\n' "$all_functions_up_to_core" | grep -qxF "$function"; then
      # Get functions after the meta ones.
      if ! printf '%s\n' "$initial_functions_plus_meta" | grep -qxF "$function"; then
        printf '%s\n' "$function"
      fi
    fi
  done
}

function list_meta_tests {
  # Want to list these tests in this order (as long as they exist).
  for test in all active inactive core; do
    if list_functions | grep -qxF "$test"; then
      printf '%s\n' "$test"
    fi
  done
  # Then programmatically list the rest.
  list_functions | while read function; do
    # Select only actual tests.
    if ! printf '%s\n' "$initial_functions" | grep -qxF "$function"; then
      # Select only meta tests.
      if printf '%s\n' "$initial_functions_plus_meta" | grep -qxF "$function"; then
        # But filter out the fixed-order ones we already listed.
        if ! printf 'all\nactive\ninactive\ncore\n' | grep -qxF "$function"; then
          printf '%s\n' "$function"
        fi
      fi
    fi
  done
}

function unit {
  unit_tests
}

function unit_tests {
  echo -e "\tUnit Tests"
  "$dirname/unit-tests.py"
}

# Capture a list of all functions defined before the tests, to tell which are actual functions
# and which are tests.
initial_functions=$(list_functions)


########## Meta tests ##########

# Run all tests.
function all {
  active
  inactive
}

function active {
  for test in $(list_active_tests); do
    $test
  done
  unit
}

function inactive {
  for test in $(list_inactive_tests); do
    $test
  done
}

function core {
  for test in $(list_core_tests); do
    $test
  done
}

# Run the errstats.py-specific tests.
function errstats {
  errstats_simple
  errstats_indels
  errstats_overlap
}

function varylen {
  varylen_barcodes
  varylen_align
  varylen_consensi
}

# Run the make-consensi.py-specific tests.
function consensi_all {
  declare -a tests
  i=1
  while read declare f test; do
    if echo "$test" | grep -qE '^consensi' && [[ $test != consensi_all ]]; then
      tests[$i]=$test
      i=$((i+1))
    fi
  done < <(declare -F)
  for test in ${tests[@]}; do
    $test
  done
}

# Get the list of functions now that the meta tests have been declared.
initial_functions_plus_meta=$(list_functions)


########## Functional tests ##########

# make-barcodes.awk
function barcodes {
  echo -e "\t${FUNCNAME[0]}:\tmake-barcodes.awk ::: families.raw_[12].fq"
  if ! local_prefix=$(_get_local_prefix "$cmd_prefix" make-barcodes.awk); then return 1; fi
  paste "$dirname/families.raw_1.fq" "$dirname/families.raw_2.fq" \
    | paste - - - - \
    | awk -f "${local_prefix}make-barcodes.awk" -v TAG_LEN=12 -v INVARIANT=5 \
    | sort \
    | diff -s - "$dirname/families.sort.tsv"
}

# make-families.sh smoke test
function make_smoke {
  if ! local_prefix=$(_get_local_prefix "$cmd_prefix" make-barcodes.awk); then return 1; fi
  echo -e "\t${FUNCNAME[0]}:\tmake-families.sh ::: smoke_[12].fq"
  "${local_prefix}make-families.sh" -t 4 -i 1 "$dirname/smoke_1.fq" "$dirname/smoke_2.fq" \
    | diff -s - "$dirname/smoke.families.tsv"
  "${local_prefix}make-families.sh" -t 4 -i 0 "$dirname/smoke_1.fq" "$dirname/smoke_2.fq" \
    | diff -s - "$dirname/smoke.families.i0.tsv"
}

# align-families.py
function align {
  echo -e "\t${FUNCNAME[0]}:\talign-families.py ::: families.sort.tsv:"
  if ! local_prefix=$(_get_local_prefix "$cmd_prefix" align-families.py); then return 1; fi
  "${local_prefix}align-families.py" --no-check-ids -q "$dirname/families.sort.tsv" \
    | diff -s - "$dirname/families.msa.tsv"
}

# align-families.py with 3 processes
function align_p3 {
  echo -e "\t${FUNCNAME[0]}:\talign-families.py -p 3 ::: families.sort.tsv:"
  if ! local_prefix=$(_get_local_prefix "$cmd_prefix" align-families.py); then return 1; fi
  "${local_prefix}align-families.py" --no-check-ids -q -p 3 "$dirname/families.sort.tsv" \
    | diff -s - "$dirname/families.msa.tsv"
}

# align-families.py smoke test
function align_smoke {
  echo -e "\t${FUNCNAME[0]}:\talign-families.py ::: smoke.families.tsv:"
  if ! local_prefix=$(_get_local_prefix "$cmd_prefix" align-families.py); then return 1; fi
  "${local_prefix}align-families.py" -q "$dirname/smoke.families.tsv" \
    | diff -s - "$dirname/smoke.families.aligned.tsv"
}

# make-consensi.py defaults on toy data
function consensi {
  _consensi families.msa.tsv families.sscs_1.fa families.sscs_2.fa families.dcs_1.fa \
            families.dcs_2.fa
}

# make-consensi.py with 3 processes
function consensi_p3 {
  _consensi families.msa.tsv families.sscs_1.fa families.sscs_2.fa families.dcs_1.fa \
            families.dcs_2.fa --processes 3
}

# make-consensi.py quality score consideration
function consensi_qual {
  _consensi qual.msa.tsv qual.10.sscs_1.fa qual.10.sscs_2.fa empty.txt empty.txt --qual 10
  _consensi qual.msa.tsv qual.20.sscs_1.fa qual.20.sscs_2.fa empty.txt empty.txt --qual 20
}

function consensi_gapqual {
  _consensi gapqual.msa.tsv gapqual.sscs_1.fa gapqual.sscs_2.fa empty.txt empty.txt --qual 25
}

function consensi_consthres {
  _consensi cons.thres.msa.tsv cons.thres.0.5.sscs_1.fa cons.thres.0.5.sscs_2.fa \
            cons.thres.0.5.dcs_1.fa cons.thres.0.5.dcs_2.fa \
            --min-cons-reads 3 --cons-thres 0.5
  _consensi cons.thres.msa.tsv cons.thres.0.7.sscs_1.fa cons.thres.0.7.sscs_2.fa \
            cons.thres.0.7.dcs_1.fa cons.thres.0.7.dcs_2.fa \
            --min-cons-reads 3 --cons-thres 0.7
}

function consensi_thres {
  _consensi cons.thres2.msa.tsv cons.thres2.sscs_1.fa cons.thres2.sscs_2.fa empty.txt empty.txt
}

function consensi_regression {
  _consensi regression.msa.tsv regression.sscs_1.fq regression.sscs_2.fq empty.txt empty.txt \
            --fastq-out 40
}

# variable-length reads
# make-barcodes.awk
function varylen_barcodes {
  echo -e "\t${FUNCNAME[0]}:\tmake-barcodes.awk ::: varylen.raw_[12].fq"
  if ! local_prefix=$(_get_local_prefix "$cmd_prefix" make-barcodes.awk); then return 1; fi
  paste "$dirname/varylen.raw_1.fq" "$dirname/varylen.raw_2.fq" \
    | paste - - - - \
    | awk -f "${local_prefix}make-barcodes.awk" -v TAG_LEN=12 -v INVARIANT=5 \
    | sort \
    | diff -s - "$dirname/varylen.sort.tsv"
}

# align-families.py
function varylen_align {
  echo -e "\t${FUNCNAME[0]}:\talign-families.py ::: varylen.sort.tsv:"
  if ! local_prefix=$(_get_local_prefix "$cmd_prefix" align-families.py); then return 1; fi
  "${local_prefix}align-families.py" --no-check-ids -q "$dirname/varylen.sort.tsv" \
    | diff -s - "$dirname/varylen.msa.tsv"
}

# make-consensi.py
function varylen_consensi {
  _consensi varylen.msa.tsv varylen.sscs_1.fa varylen.sscs_2.fa varylen.dcs_1.fa varylen.dcs_2.fa
}

# baralign.sh
function baralign {
  echo -e "\t${FUNCNAME[0]}:\tbaralign.sh ::: correct.families.d1.tsv:"
  if ! local_prefix=$(_get_local_prefix "$cmd_prefix" baralign.sh); then return 1; fi
  "${local_prefix}baralign.sh" "$dirname/correct.families.d1.tsv" "$dirname/refdir.tmp" 2>/dev/null \
    | _clean_sam | diff -s - "$dirname/correct.d1.sam"
  rm -rf "$dirname/refdir.tmp"
}

# correct.py
function correct {
  echo -e "\t${FUNCNAME[0]}:\tcorrect.py ::: correct.sam"
  if ! local_prefix=$(_get_local_prefix "$cmd_prefix" correct.py); then return 1; fi
  "${local_prefix}correct.py" --no-check-ids "$dirname/correct.families.tsv" \
      "$dirname/correct.barcodes.fa" "$dirname/correct.sam" \
    | diff -s "$dirname/correct.families.corrected.tsv" -
}

# correct.py --dist 1
function correct_d1 {
  echo -e "\t${FUNCNAME[0]}:\tcorrect.py -d 1 ::: correct.sam"
  if ! local_prefix=$(_get_local_prefix "$cmd_prefix" correct.py); then return 1; fi
  "${local_prefix}correct.py" --no-check-ids --dist 1 "$dirname/correct.families.d1.tsv" \
      "$dirname/correct.barcodes.d1.fa" "$dirname/correct.d1.sam" \
    | diff -s "$dirname/correct.families.corrected.d1.tsv" -
}


# Only the tests above are for the "core" pipeline scripts.
all_functions_up_to_core=$(list_functions)


function precheck {
  echo -e "\t${FUNCNAME[0]}:\tprecheck.py ::: families.raw_[12].fq"
  if ! local_prefix=$(_get_local_prefix "$cmd_prefix" utils/precheck.py); then return 1; fi
  "${local_prefix}precheck.py" "$dirname/families.raw_1.fq" "$dirname/families.raw_2.fq" \
    | diff -s - "$dirname/families.precheck.tsv"
}

function stats_diffs {
  echo -e "\t${FUNCNAME[0]}:\tstats.py diffs ::: gaps.msa.tsv:"
  if ! local_prefix=$(_get_local_prefix "$cmd_prefix" utils/stats.py); then return 1; fi
  "${local_prefix}stats.py" diffs "$dirname/gaps.msa.tsv" \
    | diff -s - "$dirname/gaps-diffs.out.tsv"
}


function dunovo {
  echo -e "\t${FUNCNAME[0]}:\tdunovo.py ::: families.raw_[12].fq"
  mkdir "$dirname/dunovo.tmp"
  if ! local_prefix=$(_get_local_prefix "$cmd_prefix" dunovo.py); then return 1; fi
  "${local_prefix}dunovo.py" "$dirname/families.raw_1.fq" "$dirname/families.raw_2.fq" -I \
    --min-length 20 -l "$dirname/dunovo.tmp/logs" -o "$dirname/dunovo.tmp"
  diff -s "$dirname/families.dunovo.duplex_1.fq" "$dirname/dunovo.tmp/duplex.filt_1.fq"
  diff -s "$dirname/families.dunovo.duplex_2.fq" "$dirname/dunovo.tmp/duplex.filt_2.fq"
  rm -rf "$dirname/dunovo.tmp"
}


# filter_barcodes.py
function filt {
  echo -e "\t${FUNCNAME[0]}:\tfilter_barcodes.py ::: filter.in.tsv"
  if ! local_prefix=$(_get_local_prefix "$cmd_prefix" utils/precheck.py); then return 1; fi
  "${local_prefix}filter_barcodes.py" --invalid-bases N --repeats 3 "$dirname/filter.in.tsv" \
    | diff -s "$dirname/filter.-iN-r3.out.tsv" -
  "${local_prefix}filter_barcodes.py" --invalid-bases N "$dirname/filter.in.tsv" \
    | diff -s "$dirname/filter.-iN.out.tsv" -
  "${local_prefix}filter_barcodes.py" --repeats 5 "$dirname/filter.in.tsv" \
    | diff -s "$dirname/filter.-r5.out.tsv" -
}


function errstats_simple {
  echo -e "\t${FUNCNAME[0]}:\terrstats.py ::: families.msa.tsv:"
  if ! local_prefix=$(_get_local_prefix "$cmd_prefix" utils/errstats.py); then return 1; fi
  "${local_prefix}errstats.py" --mate1 "$dirname/families.msa.tsv" \
    | diff -s "$dirname/errstats.out.tsv" -
  "${local_prefix}errstats.py" --mate1 --out-format errors1 "$dirname/families.msa.tsv" \
    | diff -s "$dirname/errstats.-R.out.tsv" -
  "${local_prefix}errstats.py" --alignment "$dirname/families.msa.tsv" \
    | diff -s "$dirname/errstats.-a.out.tsv" -
  "${local_prefix}errstats.py" --var-columns errors --alignment "$dirname/families.msa.tsv" \
    | diff -s "$dirname/errstats.-R.-a.out.tsv" -
  "${local_prefix}errstats.py" --mate1 --out-format errors2 "$dirname/families.msa.tsv" \
    | diff -s "$dirname/errstats.errors2.out.tsv" -
}

function errstats_indels {
  echo -e "\t${FUNCNAME[0]}:\terrstats.py ::: families.unequal.msa.tsv:"
  if ! local_prefix=$(_get_local_prefix "$cmd_prefix" utils/errstats.py); then return 1; fi
  < "$dirname/families.unequal.msa.tsv" "${local_prefix}errstats.py" --out-format reads --mate1 \
    | diff -s "$dirname/errstats.indels.reads.out.tsv" -
  < "$dirname/families.unequal.msa.tsv" "${local_prefix}errstats.py" --out-format reads --mate1 \
      --no-indels \
    | diff -s "$dirname/errstats.indels.reads.-I.out.tsv" -
  < "$dirname/families.unequal.msa.tsv" "${local_prefix}errstats.py" --out-format errors1 --mate1 \
    | diff -s "$dirname/errstats.indels.errors1.out.tsv" -
  < "$dirname/families.unequal.msa.tsv" "${local_prefix}errstats.py" --out-format errors1 --mate1 \
      --no-indels \
    | diff -s "$dirname/errstats.indels.errors1.-I.out.tsv" -
}


# All tests below here are considered inactive.
all_functions_minus_inactive=$(list_functions)

function errstats_overlap {
  # Note: Currently, the correct PYBAMPATH is ~/bx/code/indels/pyBamParser/hg/lib.
  # Test input was created with these commands:
  # $ ../bfx/parse-test-align.py --duplex overlap.align.txt --ref overlap.ref.fa \
  #   --fq1 overlap.reads_1.fq --fq2 overlap.reads_2.fq
  # $ ../bfx/align-bowtie.sh -c -b overlap.reads.bam overlap.ref.fa overlap.reads_[12].fq
  echo -e "\t${FUNCNAME[0]}:\terrstats.py ::: overlap.families.msa.tsv"
  if ! local_prefix=$(_get_local_prefix "$cmd_prefix" utils/errstats.py); then return 1; fi
  "${local_prefix}errstats.py" --dedup --min-reads 3 --bam "$dirname/overlap.sscs.bam" \
    "$dirname/overlap.families.msa.tsv" --overlap-stats "$dirname/overlap.overlaps.tmp.tsv" >/dev/null
  diff -s "$dirname/overlap.overlaps.tmp.tsv" "$dirname/overlap.overlaps.expected.tsv"
  if [[ -f "$dirname/overlap.overlaps.tmp.tsv" ]]; then
    rm "$dirname/overlap.overlaps.tmp.tsv"
  fi
}

# utility function for all make-consensi.py tests
function _consensi {
  # Read required arguments.
  input=$1
  sscs1=$2
  sscs2=$3
  dcs1=$4
  dcs2=$5
  # Read optional arguments (after the required ones).
  declare -a args
  i=6
  while [[ ${!i} ]]; do
    args[$i]=${!i}
    i=$((i+1))
  done
  echo -e "\t${FUNCNAME[1]}:\tmake-consensi.py ${args[@]} ::: $input:"
  if ! local_prefix=$(_get_local_prefix "$cmd_prefix" make-consensi.py); then return 1; fi
  "${local_prefix}make-consensi.py" --aligner biopython ${args[@]} "$dirname/$input" \
    --sscs1 "$dirname/cons.tmp.sscs_1.fa" --sscs2 "$dirname/cons.tmp.sscs_2.fa" \
    --dcs1  "$dirname/cons.tmp.dcs_1.fa"  --dcs2  "$dirname/cons.tmp.dcs_2.fa"
  diff -s "$dirname/cons.tmp.sscs_1.fa" "$dirname/$sscs1"
  diff -s "$dirname/cons.tmp.sscs_2.fa" "$dirname/$sscs2"
  diff -s "$dirname/cons.tmp.dcs_1.fa"  "$dirname/$dcs1"
  diff -s "$dirname/cons.tmp.dcs_2.fa"  "$dirname/$dcs2"
  for file in cons.tmp.sscs_1.fa cons.tmp.sscs_2.fa cons.tmp.dcs_1.fa cons.tmp.dcs_2.fa; do
    if [[ -f "$dirname/$file" ]]; then
      rm "$dirname/$file"
    fi
  done
}

function _clean_sam {
  # Remove @PG line and XM:i: tags.
  awk -F '\t' -v OFS='\t' '
    $1 !~ /^@PG$/ {
      for (i=1; i<=NF; i++) {
        if (i == 1) {
          printf("%s", $i)
        } else if (i <= 11 || substr($i, 1, 5) != "XM:i:") {
          printf("\t%s", $i)
        }
      }
      printf("\n")
    }'
}

function _get_local_prefix {
  local cmd_prefix="$1"
  postfix="$2"
  base=$(basename "$postfix")
  if [[ "$cmd_prefix" ]]; then
    path="${cmd_prefix}$postfix"
    local_prefix=$(dirname "${cmd_prefix}$postfix")/
  else
    # If $cmd_prefix is blank, the user wants to try to execute scripts via the $PATH.
    # Use the basename of the $postfix, removing any directories before the actual script.
    path="$base"
    local_prefix=
  fi
  if which "$path" >/dev/null 2>/dev/null; then
    echo "$local_prefix"
  elif [[ -f "$path" ]]; then
    echo "$local_prefix"
  else
    echo -e "\e[31mError: $base missing!\e[m Searched for: \"$path\"" >&2
    return 1
  fi
}

main "$@"
