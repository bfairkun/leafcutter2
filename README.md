# LeafCutter2

## Introduction

LeafCutter2 is a tool for clustering, functional characterization, and quantification of splice junction read counts. It implements a novel dynamic programming algorithm over the standard LeafCutter output to classify splice junctions and alternative splicing events according to their molecular function (i.e.: protein-coding or unproductive).

[LeafCutter2 also provides a Python implementation of our original differential analysis algorithm](https://github.com/leafcutter2/leafcutter-ds). This is faster and easier to install than the R version, and it is fully compatible with our new clustering and functional classification tools.

## Installation

Because the example data is bundled with the repository, we recommend cloning from GitHub and installing from the local clone:

```bash
git clone https://github.com/leafcutter2/leafcutter2.git
cd leafcutter2
```

All Python dependencies can be installed into a conda environment using the bundled environment file:

```bash
conda env create -f leafcutter2_env.yml
conda activate leafcutter2
```

The environment can also be created using [Mamba](https://mamba.readthedocs.io/en/latest/), which is generally faster:

```bash
mamba env create -f leafcutter2_env.yml
conda activate leafcutter2
```

Then install LeafCutter2 itself:

```bash
pip install .
```

After installation the following commands will be available on your PATH: `leafcutter2`, `leafcutter2-make-clusters`, and `leafcutter2-star2junc`.


## Clustering, classifying and quantifying splice junctions

### Input:
- Splice junction BED files. They should contain at least six columns (chrom, start, end, name, score and strand), with the fifth column corresponding to the splice junction read counts.
- GTF annotation with genes, start codons and stop codons.
- Genome assembly FASTA file. It must correspond to the same assembly as the GTF file. Make sure it has an index in a [faidx file](https://www.htslib.org/doc/samtools-faidx.html) in the same directory.

We recommend using the BED-formatted `.junc` files obtained from BAM files using [regtools junctions extract](https://regtools.readthedocs.io/en/latest/commands/junctions-extract/) as input.

E.g.: `regtools junctions extract -a 8 -i 50 -I 500000 bamfile.bam -o outfile.junc`.

The BED files can also be obtained from [STAR's](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) `SJ.out.tab` files using the bundled conversion tool:

```bash
leafcutter2-star2junc sample.SJ.out.tab sample.bed.gz
```

### Running LeafCutter2

We provide a basic example from GTEx in the `example/` directory. From within that directory, a basic LeafCutter2 run works as follows:

```bash
cd example/
leafcutter2 \
    -j junction_files.txt \
    -r output_dir \
    -A annotation/chr10.gtf.gz \
    -G annotation/chr10.fa.gz
```

-    `-j junction_files.txt` should be a text file listing paths to each junction file, one path per line.
-    `-r output_dir` specifies the directory of output (default is current directory, `./`).
-    `-A annotation/chr10.gtf.gz` is a gtf file of chromosome 10 obtained from Gencode v43.
-    `-G annotation/chr10.fa.gz` a FASTA file of chromosome 10 (GRCh38 assembly).

**Note:** Make sure that the chromosome names match between the BED, GTF and FASTA files.

### Output:
- `leafcutter2.cluster_ratios.gz` a table quantifying splice junction read counts for each intron, divided by the total number of reads observed for the intron cluster to which the intron belongs. Each row corresponds to a splice junction, with the first column indicating the splice junction ID and all subsequent columns corresponding to each sample. The splice junction ID has the following format: `chr10:134786:179993:clu_1_+:PR`, indicating the chromosome, start and end of the splice junction, the LeafCutter intron cluster to which it belongs (in this case, `clu_1_+`), and a label indicating the splice junction's function: **PR** (productive/protein-coding), **UP** (unproductive), **NE** (ambiguous in their functional effect) or **IN** (intergenic).
- `leafcutter2.junction_counts.gz` a table similar to `leafcutter2.cluster_ratios.gz`, except it only keeps track of the splice junction read counts for each intron. Useful for differential splicing analysis with [LeafCutter's R package](https://davidaknowles.github.io/leafcutter/).
- `clustering/` a directory containing files relevant for clustering and annotation.
    - `clustering/leafcutter2_clusters` contains the intron clusters. Useful for skipping the clustering step in repeated runs.
    - Other files documenting stats from the clustering and classification algorithm. Useful for debugging.

### A note on GTF files:

For classification, LeafCutter2 requires information from a GTF file. Please ensure that the features column (third column on a GTF) includes the following type of information, and that the names match:
- `gene`
- `transcript`
- `CDS`
- `start_codon`
- `stop_codon`

In addition, LeafCutter2 requires information on gene and transcript types. These come in the 9th column of a GTF file, but different annotations can use different tags (e.g., Gencode typically uses `transcript_type`, while Ensembl uses `transcript_biotype`). You can specify what tag your GTF uses by using the following parameters:
- `--transcript_type` (default: auto-detect)

You can also specify what tag your GTF uses for gene and transcript name (or if you prefer, change from names to gene and transcript IDs) by using:
- `--gene_name` (default: `gene_name`)
- `--transcript_name` (default: `transcript_name`)

Finally, some GTFs might lack some or all of these features and information (e.g., assembled GTFs from StringTie, or some GTFs from the UCSC Genome Browser). LeafCutter2 will attempt to reformat such GTFs automatically. You can disable this with `--no-auto-reformat-gtf`.


### Parameters

```
leafcutter2 --help

usage: leafcutter2 [-h] (-j JUNCFILES | --leafcutter1-counts-file COUNTS_FILE)
                   [-o OUTPREFIX] [-q] [-r RUNDIR] [-l MAXINTRONLEN]
                   [-m MINCLUREADS] [-M MINREADS] [-D MINREADSTD]
                   [-p MINCLURATIO] [-c CLUSTER] [-k] [-C] [-A ANNOT]
                   [-G GENOME] [-f OFFSET] [-T] [-L] [-P] [-t TRANSCRIPT_TYPE]
                   [-gn GENE_NAME] [-tn TRANSCRIPT_NAME] [-N N]
                   [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
                   [--no-auto-reformat-gtf]

options:
  -h, --help            show this help message and exit
  -j JUNCFILES, --juncfiles JUNCFILES
                        a text file storing paths to junction files, one path per line
  --leafcutter1-counts-file COUNTS_FILE
                        pre-existing counts file (*.counts.gz). If provided, skip
                        clustering and counting steps and proceed directly to classification
  -o OUTPREFIX, --outprefix OUTPREFIX
                        output prefix (default leafcutter2)
  -q, --quiet           don't print status messages to stdout
  -r RUNDIR, --rundir RUNDIR
                        write to directory (default ./)
  -l MAXINTRONLEN, --maxintronlen MAXINTRONLEN
                        maximum intron length in bp (default 100,000bp)
  -m MINCLUREADS, --minclureads MINCLUREADS
                        minimum reads in a cluster (default 30 reads)
  -M MINREADS, --minreads MINREADS
                        minimum reads for a junction to be considered for
                        clustering (default 5 reads)
  -D MINREADSTD, --minreadstd MINREADSTD
                        minimum standard deviation of reads across samples for a
                        junction to be included in output (default 0.5)
  -p MINCLURATIO, --mincluratio MINCLURATIO
                        minimum fraction of reads in a cluster that support a
                        junction (default 0.001)
  -c CLUSTER, --cluster CLUSTER
                        skip intron clustering step, use pre-determined clusters
  -k, --checkchrom      check that the chromosomes are well formatted e.g. chr1,
                        chr2, ..., or 1, 2, ...
  -C, --includeconst    also include constitutive introns
  -A ANNOT, --annot ANNOT
                        GTF annotation file
  -G GENOME, --genome GENOME
                        Genome fasta file
  -N N, --max_juncs N   skip solveNMD function if gene contains more than N juncs.
                        Juncs in skipped genes are assigned Coding=False. Default 10000
  -f OFFSET, --offset OFFSET
                        Offset sometimes useful for off by 1 annotations. (default 0)
  -T, --keeptemp        keep temporary files. (default false)
  -L, --keepleafcutter1
                        keep LeafCutter1-compatible files. Useful for running
                        differential splicing analysis with leafcutter's R package.
                        (default false)
  -P, --keepannot       save parsed annotations to .pckle files. (default false)
  -t TRANSCRIPT_TYPE, --transcript_type TRANSCRIPT_TYPE
                        tag for transcript type in GTF file (default: auto-detect
                        from transcript_type, transcript_biotype)
  -gn GENE_NAME, --gene_name GENE_NAME
                        tag for gene name or ID in GTF file (default: auto-detect
                        from gene_name, gene_id, gene_symbol)
  -tn TRANSCRIPT_NAME, --transcript_name TRANSCRIPT_NAME
                        tag for transcript name or ID in GTF file (default: auto-detect
                        from transcript_name, transcript_id)
  --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        logging level (default INFO)
  --no-auto-reformat-gtf
                        disable automatic GTF reformatting when validation fails
```

## Pre-clustering splice junctions (optional)

Generating intron clusters first can save time for multiple subsequent runs. The `leafcutter2-make-clusters` command makes intron clusters separately that can be later used as input for `leafcutter2`. You can generate clusters by running:

```bash
leafcutter2-make-clusters \
    -j junction_files.txt \
    -r output_dir
```

This will generate a file named `output_dir/clustering/leafcutter2_clusters` that can be later used as an input for LeafCutter2 to skip the clustering step:

```bash
leafcutter2 \
    -j junction_files.txt \
    -r output_dir \
    -A annotation/chr10.gtf.gz \
    -G annotation/chr10.fa.gz \
    -c output_dir/clustering/leafcutter2_clusters
```

**Note:** The junction-filtering options (`--minclureads`, `--minreads`, `--mincluratio`) will be ignored by `leafcutter2` if a pre-defined set of intron clusters is provided.
