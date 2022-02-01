#!/usr/bin/env bash
if [ x$BASH = x ] || [ ! $BASH_VERSINFO ] || [ $BASH_VERSINFO -lt 4 ]; then
  echo "Error: Must use bash version 4+." >&2
  exit 1
fi
set -ue

Project=dunovo
DefaultChunkMbs=512
RefdirDefault=refdir
EtDomain='nstoler.com'
RequiredCommands='bowtie bowtie-build samtools awk'
verbosity=2

Usage="Usage: \$ $(basename $0) [options] families.tsv [refdir [outfile.sam|outfile.bam]]
families.tsv: The families file produced by make-barcodes.awk and sorted.
refdir:  The directory to put the reference file (\"barcodes.fa\") and its index
         files in. Default: \"$RefdirDefault\".
outfile: Print the output to this path. It will be in SAM format unless the
         path ends in \".bam\". If not given, it will be printed to stdout
         in SAM format.
-R: Don't include reversed barcodes (alpha+beta -> beta+alpha) in the alignment
    target.
-t: Number of threads for bowtie and bowtie-build to use (default: 1).
-c: Number to pass to bowtie's --chunkmbs option (default: $DefaultChunkMbs).
-p: Report helpful usage data to the developer, to better understand the use
    cases and performance of the tool. The only data which will be recorded is
    the name and version of the tool, the size of the input data, the time taken
    to process it, the IP address of the machine running it, and some
    performance-related parameters (-t, -c, and the format of the output file).
    No filenames are sent. All the reporting and recording code is available at
    https://github.com/NickSto/ET.
-g: Report the platform as \"galaxy\" when sending usage data.
-v: Print the version number and exit.
-q: Quiet mode"

function main {

  script_dir=$(get_script_dir)
  version=$(version "$script_dir")
  start_time=$(date +%s)

  # Read in arguments and check them.
  if [[ "$#" -ge 1 ]] && [[ "$1" == '--version' ]]; then
    echo "$version"
    return
  fi

  threads=1
  reverse=true
  chunkmbs="$DefaultChunkMbs"
  phone=
  platform_args=
  while getopts "rhc:t:pgvq" opt; do
    case "$opt" in
      r) reverse='';;
      t) threads="$OPTARG";;
      c) chunkmbs="$OPTARG";;
      p) phone='home';;
      g) platform_args='--platform galaxy';;
      v) echo "$version" && return;;
      q) verbosity=1;;
      [h?]) fail "$Usage";;
    esac
  done
  # Get positional arguments.
  families="${@:$OPTIND:1}"
  refdir="${@:$OPTIND+1:1}"
  outfile="${@:$OPTIND+2:1}"

  if [[ "$phone" ]] && [[ -x "$script_dir/ET/phone.py" ]]; then
    #TODO: Use version.py --get-key to read the project from VERSION.
    set +e
    run_id=$("$script_dir/ET/phone.py" start --domain "$EtDomain" \
             --project "$Project" --script "$(basename "$0")" \
             --version "$version" $platform_args)
    set -e
  fi

  # Validate arguments.
  if ! [[ "$families" ]]; then
    fail "$Usage"$'\n\nError: Must provide an input families.tsv file.'
  elif ! [[ -f "$families" ]]; then
    fail "Error: families_file \"$families\" not found."
  fi
  if ! [[ "$refdir" ]]; then
    refdir="$RefdirDefault"
  fi
  if ! [[ -d "$refdir" ]]; then
    echo "Info: ref_dir \"$refdir\" not found. Creating.." >&2
    mkdir "$refdir"
  fi
  # Determine how and where to put the output.
  if [[ "${outfile:${#outfile}-4}" == .bam ]]; then
    format=bam
  else
    format=sam
  fi
  sam_outfile=
  outbase=$(echo "$outfile" | sed -E 's/\.bam$//')
  if [[ "$outfile" ]]; then
    if [[ -e "$outfile" ]]; then
      fail "Error: output file \"$outfile\" already exists."
    fi
    if [[ "$format" == bam ]]; then
      if [[ -e "$outbase.sam" ]] || [[ -e "$outbase.bam.bai" ]]; then
        fail "Error: A related filename already exists (.sam/.bam.bai)."
      fi
      sam_outfile="$outbase.sam"
    else
      sam_outfile="$outfile"
    fi
  fi

  # Check for required commands.
  for cmd in $RequiredCommands; do
    if ! which "$cmd" >/dev/null 2>/dev/null; then
      fail "Error: command \"$cmd\" not found."
    fi
  done

  # Check version of bowtie-build.
  # Only version 1.2.1 and above had --threads option.
  indexer_is_threaded=$(bowtie-build --version | awk '
    $1 == "bowtie-build" && $2 == "version" {
      split($3, fields, ".")
      maj_min = fields[1] "." fields[2]
      if (maj_min > 1.2) {
        print "yes"
      } else if (maj_min == 1.2 && fields[3] >= 1) {
        print "yes"
      }
    }')
  if [[ "$indexer_is_threaded" ]]; then
    indexer_threads="--threads $threads"
  else
    indexer_threads=
  fi

  if [[ "$phone" ]] && [[ -x "$script_dir/ET/phone.py" ]]; then
    set +e
    size=$(du -sb "$families" | awk '{print $1}')
    run_data="\"format\":\"$format\", \"threads\":\"$threads\", \"chunkmbs\":\"$chunkmbs\",\
              \"families_size\":\"$size\""
    "$script_dir/ET/phone.py" prelim --domain "$EtDomain" \
      --project "$Project" --script "$(basename "$0")" \
      --version "$version" $platform_args --run-id "$run_id" \
      --run-data "{$run_data}"
    set -e
  fi

  if [[ "$verbosity" -ge 2 ]]; then
    printf 'families: %s\n' "$families" >&2
    printf 'refdir:   %s\n' "$refdir" >&2
    printf 'format:   %s\n' "$format" >&2
    printf 'outfile:  %s\n' "$outfile" >&2
    printf 'outbase:  %s\n' "$outbase" >&2
  fi

  # Create FASTA with barcodes as "reads" for alignment.
  awk '$1 != last {
    count++
    print ">" count
    print $1
  }
  {
    last = $1
  }' "$families" > "$refdir/barcodes.fa"

  # Create "reference" to align the barcodes to.
  if [[ "$reverse" ]]; then
    # If we're including reversed barcodes, create a new FASTA which includes reversed barcodes
    # as well as their forward versions.
    awk '
      $1 != last {
        count++
        bar = $1
        print ">" count
        print bar
        print ">" count ":rev"
        print swap_halves(bar)
      }
      {
        last = $1
      }
      function swap_halves(str) {
        half = length(str)/2
        alpha = substr(str, 1, half)
        beta = substr(str, half+1)
        return beta alpha
      }' "$families" > "$refdir/barcodes-ref.fa"
  else
    # If we're not including reversed barcodes, the original FASTA is all we need. Just link to it.
    ln -s "$refdir/barcodes.fa" "$refdir/barcodes-ref.fa"
  fi

  # Perform alignment.
  exho bowtie-build -f $indexer_threads --offrate 1 \
    "$refdir/barcodes-ref.fa" "$refdir/barcodes-ref" > "$refdir/bowtie-build.out.log"
  set +e
  exho bowtie --chunkmbs "$chunkmbs" --threads "$threads" -f --sam -a --best -v 3 \
    "$refdir/barcodes-ref" "$refdir/barcodes.fa" "$sam_outfile" 2> "$refdir/bowtie.err.log"
  exit_code="$?"
  set -e
  if [[ "$verbosity" -ge 2 ]]; then
    cat "$refdir/bowtie.err.log" >&2
  fi
  if [[ "$exit_code" != 0 ]]; then
    exit "$exit_code"
  fi
  if [[ "$outfile" ]] && [[ "$format" == bam ]]; then
    exho samtools view -Sb "$sam_outfile" | exho samtools sort -o - dummy > "$outfile"
    if [[ -s "$outfile" ]]; then
      exho samtools index "$outfile"
      exho rm "$sam_outfile"
    fi
  fi
  # Check output.
  success=null
  if [[ "$outfile" ]]; then
    if [[ -s "$outfile" ]]; then
      if [[ "$format" == bam ]] && [[ -e "$outbase.sam" ]]; then
        exho rm "$outbase.sam"
      fi
      success=true
      if [[ "$verbosity" -ge 2 ]]; then
        echo "Success. Output located in \"$outfile\"." >&2
      fi
    else
      success=false
      fail "Warning: No output file \"$outfile\" found."
    fi
  fi

  if [[ "$phone" ]] && [[ -x "$script_dir/ET/phone.py" ]]; then
    set +e
    now=$(date +%s)
    run_time=$((now-start_time))
    "$script_dir/ET/phone.py" end --domain "$EtDomain" \
      --project "$Project" --script "$(basename "$0")" \
      --version "$version" $platform_args --run-id "$run_id" --run-time "$run_time" \
      --run-data "{$run_data, \"success\":$success}"
    set -e
  fi
}

function version {
  if [[ "$#" -ge 1 ]]; then
    script_dir="$1"
  else
    script_dir=$(get_script_dir)
  fi
  if [[ -x "$script_dir/utillib/version.py" ]]; then
    "$script_dir/utillib/version.py" --config-path "$script_dir/VERSION" --repo-dir "$script_dir"
  fi
}

function get_script_dir {
  # Find the actual directory this file resides in (resolving links).
  if readlink -f dummy >/dev/null 2>/dev/null; then
    script_path=$(readlink -f "${BASH_SOURCE[0]}")
  else
    # readlink -f doesn't work on BSD systems.
    script_path=$(perl -MCwd -le 'print Cwd::abs_path(shift)' "${BASH_SOURCE[0]}")
  fi
  dirname "$script_path"
}

function exho {
  # Echo a command to stderr, then execute it.
  if [[ "$verbosity" -ge 2 ]]; then
    echo "$@" >&2
  fi
  "$@"
}

function fail {
  echo "$@" >&2
  exit 1
}

main "$@"
