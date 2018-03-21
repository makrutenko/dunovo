# Du Novo

This is a pipeline for processing duplex sequencing data without the use of a reference genome.

The pipeline was designed for use with the duplex method described in [Kennedy *et al.* 2014](https://dx.doi.org/10.1038/nprot.2014.170), but the assumptions are relatively minimal, so you should be able to apply it to variants of the protocol.

Du Novo 2.0 is released under the GPLv2 license, except for some portions governed by the MIT license. Earlier versions were released under the BSD license. See `LICENSE.txt` for details.


## Running Du Novo from Galaxy

We created a comprehensive [tutorial](https://github.com/galaxyproject/dunovo/wiki) explaining all aspects of interactive use of Du Novo from within [Galaxy](http://usegalaxy.org).


## Running Du Novo on the command line


### Dependencies

#### Required

The pipeline requires a Unix operating system and Python version 2.7. Linux is recommended. OS X and BSD may work, but are untested.

It also requires several standard Unix tools. Version numbers in parentheses are what the software was tested with, but other versions likely work. These must be available on your [`$PATH`](https://en.wikipedia.org/wiki/Search_path):  
 -  the [`gcc`](https://gcc.gnu.org/) command (4.8.4)
 -  the [`make`](https://www.gnu.org/software/make/) command (3.81)
 -  the [`bash`](https://www.gnu.org/software/bash/bash.html) command (4.0)
 -  the [`awk`](https://www.gnu.org/software/gawk/) command (4.0.1)
 -  the [`paste`](https://www.gnu.org/software/coreutils/coreutils.html) command (8.21)
 -  the [`sort`](https://www.gnu.org/software/coreutils/coreutils.html) command (8.21)


#### Optional

To use `align-families.py`'s `-a mafft` option, this must be available on your `$PATH`:  
 - the [`mafft`](http://mafft.cbrc.jp/alignment/software/) command (v7.271 or v7.123b)

To use the barcode error correction script `correct.py`, the following modules must be available from Python and the commands must be on your `$PATH`:
 - the [networkx](https://pypi.python.org/pypi/networkx) Python module (1.9, 1.10, or 1.11)
 - the [`bowtie`](http://bowtie-bio.sourceforge.net/index.shtml) command (1.2.1.1) (nothing below 1.1.2 is confirmed to work)
 - the [`bowtie-build`](http://bowtie-bio.sourceforge.net/index.shtml) command (1.2.1.1) (same)
 - the [`samtools`](http://samtools.sourceforge.net/) command (0.1.18)


### Download

#### Git

    $ git clone --recursive https://github.com/galaxyproject/dunovo.git
    $ cd dunovo
    $ git checkout master
    $ git submodule update --recursive

#### Via GitHub webpage

Click the [releases](https://github.com/galaxyproject/dunovo/releases) tab at the top of this page, and find the latest release. Download the zip file (the "Source code (zip)" link), as well as `kalign.zip`, `utillib.zip`, and `ET.zip`. Extract the first zip file, then unzip the other three into the extracted directory and name those directories `kalign`, `utillib`, and `ET`, respectively.

In the end, the organization and names of the three directories should look like this:

    dunovo
    ├─╴kalign
    ├─╴utillib
    └─╴ET

### Installation

    $ cd dunovo
    $ make

The `make` command is needed to compile the C modules and kalign, which are required. You need to be in the root source directory (where the file `Makefile` is) before running the command.


### Usage

This example shows how to go from raw duplex sequencing data to the final duplex consensus sequences.

Your raw reads should be in `reads_1.fastq` and `reads_2.fastq`. And the scripts `align-families.py`, `make-consensi.py`, `baralign.sh`, and `correct.py` should be on your `PATH`. Also, in the following command, replace `make-barcodes.awk` with the actual path to that script (included in this pipeline).

1. Sort the reads into families based on their barcodes and split the barcodes from the sequence.  
    ```bash
    $ paste reads_1.fastq reads_2.fastq \
      | paste - - - - \
      | awk -f make-barcodes.awk \
      | sort > families.tsv
    ```

2. (Optional) Correct errors in barcodes.
    ```bash
    $ baralign.sh families.tsv refdir barcodes.bam
    $ samtools view -f 256 barcodes.bam \
      | correct.py families.tsv refdir/barcodes.fa \
      | sort > families.corrected.tsv
    ```

3. Do multiple sequence alignments of the read families.  
`$ align-families.py families.tsv > families.msa.tsv`

4. Build duplex consensus sequences from the aligned families.  
`$ make-consensi.py families.msa.tsv -1 duplex_1.fa -2 duplex_2.fa`

See all options for a given command by giving it the `-h` flag.


### Details

#### 1. Sort the reads into families based on their barcodes and split the barcodes from the sequence.  

    $ paste reads_1.fastq reads_2.fastq \
      | paste - - - - \
      | awk -f make-barcodes.awk \
      | sort > families.tsv


This command pipeline will transform each pair of reads into a one-line record, split the 12bp barcodes off them, and sort by their combined barcode. The end result is a file (named `families.tsv` above) listing read pairs, grouped by barcode. See `make-barcodes.awk` for the details on the formation of the barcodes and the format.

Note: This step requires your FASTQ files to have exactly 4 lines per read (no multi-line sequences). Also, in the output, the read sequence does not include the barcode or the 5bp constant sequence after it. You can customize the length of the barcode or constant sequence by setting the awk constants `TAG_LEN` and `INVARIANT` (i.e. `awk -v TAG_LEN=10 make-barcodes.awk`).

#### 2. (Optional) Correct errors in barcodes.

    $ baralign.sh families.tsv refdir barcodes.bam
    $ samtools view -f 256 barcodes.bam \
      | correct.py families.tsv refdir/barcodes.fa \
      | sort > families.corrected.tsv

These commands takes the `families.tsv` file produced in the previous step, "corrects"\* the barcodes in it, and outputs a new version of `families.tsv` with the new barcodes. It does this by aligning all barcodes to themselves, finding pairs of barcodes which differ by only a few edits. Grouping sets of related barcodes gives groups which are likely descended from the same original barcode, differing only because of PCR and/or sequencing errors. By default, only barcodes that differ by 1 edit are allowed. You can allow greater edit distances between barcodes with the `--dist` option to `correct.py`.

\*"corrects" is in scare quotes because the algorithm isn't actually focused on finding the original barcode sequence. Its purpose is instead to find barcodes which are all descended from the same original barcode, but now have different sequences because of errors. It finds each group of related barcodes and replaces them with a single barcode, so that the following steps identify them as one family.


#### 3. Do multiple sequence alignments of the read families.  

`$ align-families.py families.tsv > families.msa.tsv`

This step aligns each family of reads, but it processes each strand separately. It can be parallelized with the `-p` option.

By default, this uses the Kalign2 multiple sequence alignment algorithm. Use `-a mafft` to select MAFFT instead. Kalign2 is reccommended, as its results are of similar accuracy and it's 7-8x faster.


#### 4. Build duplex consensus sequences from the aligned families.  

`$ make-consensi.py families.msa.tsv -1 duplex_1.fa -2 duplex_2.fa`

This calls a consensus sequence from the multiple sequence alignments of the previous step. It does this in two steps: First, single-strand consensus sequences (SSCSs) are called from the family alignments, then duplex consensus sequences are called from pairs of SSCSs.

When calling SSCSs, by default 3 reads are required to successfully create a consensus from each strand (change this with `-r`). Quality filtering is done at this step by excluding bases below a quality threshold. By default, no base with a PHRED quality less than 20 will contribute to the consensus (change this with `-q`). If no base passes the threshold or there is no majority base, `N` will be used.

The duplex consensus sequences are created by comparing the two SSCSs. For each base, if they agree, that base will be inserted. If they disagree, the IUPAC ambiguity code for the two bases will be used. Note that a disagreement between a base and a gap will result in an `N`.

The output of this step is the duplex consensus sequences in FASTA format. To also output all single-strand consensus sequences (including those which didn't produce a duplex consensus), use the `--sscs1` and `--sscs2` options.

The reads will be printed in two files, one per paired-end mate, with this naming format:  
`>{barcode} {# reads in strand 1 family}-{# reads in strand 2 family}`  
e.g.  
`>TTGCGCCAGGGCGAGGAAAATACT 8-13`


### Known bugs

Be aware that a [known bug](https://stackoverflow.com/questions/1408356/keyboard-interrupts-with-pythons-multiprocessing-pool/1408476#1408476) in Python when using the [multiprocessing](https://docs.python.org/2/library/multiprocessing.html) module makes it impossible to kill a running process with the Ctrl+C keyboard command. So if you run `align-families.py` or `make-consensi.py` in the foreground, you'll have to exit via Ctrl+Z to stop and background the job, then kill the process (e.g. with `$ kill %1`, if it's the only backgrounded job).
