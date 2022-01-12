## Tutorial 3

**Objective:** To perform a de novo assembly of the Wuhan-1 genome using two different assemblers and compare the results.
<br>

Regardless of the sequencing technology (Sanger, 454, Illumina, PacBio, Nanopore, etc.)
and the assembly algorithm (de bruijn graph, overlap layout consensus, etc.),
all de novo genome assembly programs essentially rely on stringing together overlapping
portions of reads into longer fragments called contigs (contiguous sequences).
In turn, these contigs can be joined to generate scaffolds if you have paired-end reads (forward and reverse reads). Here
we will attempt to recover the original Wuhan-1 genome, which helped to identify SARS-CoV-2 as a novel Betacoronavirus.
SPAdes and MEGAHIT are two of the most frequently employed assemblers for bacterial/viral genomes.  MEGAHIT was designed for metagenomes
but can also work on single isolate libraries. Shovill is basically a faster adaptation of SPAdes.
<br>


## Trim adapters and quality filter

Use TrimGalore (or try BBtools!) to filter your reads.


	trim_galore -q 20 --paired SRR10971381_1.fastq SRR10971381_2.fastq

> We would hope that reads submitted to the SRA have adapters trimmed, but this is often not the case.

</br>

Rename your files something less clunky.


	mv SRR10971381_1_val_1.fq SRR10971381-trimmed_1.fq
	mv SRR10971381_2_val_2.fq SRR10971381-trimmed_2.fq


<br>

## Downsample your trimmed fastq files with seqtk

	seqtk sample -s100 SRR10971381-trimmed_1.fq 5000000 > SRR10971381-sub_1.fastq
	seqtk sample -s100 SRR10971381-trimmed_2.fq 5000000 > SRR10971381-sub_2.fastq

> [Setqk](https://github.com/lh3/seqtk) is like a swiss-army knife for manipulating fasta/fastq files.
> To ensure that we randomly subsample the same reads from both R1 and R2 files (maintain the files as paired),
> we need to use the same starting seed (specified by the '-s' argument). Subsampling is necessary because our fastq files are over 8 Gb in size (>28 million reads)!
> We don't have enough memory on our VMs to run a full assembly.
>
> Since this was a metatranscriptome from a human subject, we might expect the majority of our reads to be human-derived. We
> might want to exclude these sequences as well before assembling to decrease the run-time.  In this case, because we downloaded the fastq files
> from the SRA, all human DNA/RNA should be removed already.

<br>

Copy a slightly modified, downsampled version of the original fastqs from our bucket to your VM.


	gsutil -m cp gs://wc-bms-bi-training-bucket/reference_assembly/downsampled*fastq .

> This is a further downsampled version of the original data.  Use these fastqs for all downstream processes in this section of the tutorial.

<br>


## Perform a metagenome assembly with SPAdes


Assemble the metagenomic reads from our SRA download in Tutorial 1 with [SPAdes](https://github.com/ablab/spades).


	docker run --rm -v $(pwd):/data -w /data staphb/spades spades -1 downsampled-trimmed_1.fastq -2 downsampled-trimmed_2.fastq -o wuhan_assembly_spades/ --meta &

**Note** Now that our VMs have additional RAM, we should be able to use all the reads. Below, we are running the SPAdes installed on our VM and foregoing the read error correction step to save time.

	spades --only-assembler -1 SRR10971381-trimmed_1.fq -2 SRR10971381-trimmed_2.fq -o wuhan_assembly_spades

> SPAdes provides several assembly options (including the metagenome option) and is ideally suited for smaller genomes.
> SPAdes will perform error-correction prior to assembling or you can run without the error-correction step as we have done in the second example.
>
> **Note:** The ampersand specifies performing a job in the background. This returns your command-line prompt so that you can continue working.
> Given the size of the library (8 GB!), assembly could take quite some time if we used all reads. You can see that your job is stil running by typing `top`,
> which allows you to see the status of all jobs on your VM. But because SPAdes prints a lot of information to STDOUT, it's rather inconvenient to work
> in this terminal.  You can start a new connection by selecting this option in the drop down menu that looks like a cog in the upper right corner of your terminal.

<br>

## Generate some summary statistics for your assembly

The final assembly is located in the 'scaffolds.fasta' file in the output directory you specified.
Ideally, we would like a closed genome (one chromosome).
Although short reads are typically not enough to scaffold contigs into chromosome-sized pieces, our genome is only
29903 base pairs long.

* How many contigs are in the 'scaffolds.fasta' file?
* What is the average length of your scaffolds?  Hint: the length of each scaffold is given in the definition line.
* What is the N50 of your assembly (N50 is the value whereby 50% of an assembly is contained in contigs/scaffolds equal to or larger than this value)?
* What is the average coverage of your scaffolds? Hint: the average coverage of each scaffold is also given in the definition line.
* Since there are several contigs in this file, how can we identify which are SARS-CoV-2 (Hint: try BLAST)?
* If there are several SARS-CoV-2 contigs, how can we identify the Wuhan-1 genome?


<br>


## BLAST (Basic Local Alignment Search Tool)

[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) looks for regions of similarity between sequences, which serves as a proxy of homology (descended from a common ancestor).
We could use the website to perform a BLAST search and identify similar sequences to our assembled contigs our we could use a locally installed version.  Similarly, we could use
the nucleotide database provided by NCBI (all nucleotide sequences deposited at NCBI) or we could make our own database. Since we are only trying to identify the Wuhan-1 SARS-CoV-2 genome
from our assembly, we could simply make a database of the Wuhan-1 reference genome and BLAST to it.
Make a BLAST database of the Wuhan-1 genome downloaded from NCBI in Tutorial 1


	makeblastdb -in wuhan.fna -dbtype nucl

> Here we are specifying the input file as well as the database type, which is nucleotide.  We could also generate a database from
protein sequences, which would take the argument `dbtype prot`.

<br>

## Blast your assembly to the Wuhan-1 genome.


	blastn -query wuhan_assembly_spades/scaffolds.fasta -db wuhan.fna -outfmt 6 -out spades.br

> There are different flavors of blast depending on the query and database type (nucleotide vs. protein).
> The `blastn` command specifies a nucleotide to nucleotide search.  We can modify the format of
> the output with the `-outfmt` option.  Here we are saving our output to a tab-delimited file called 'spades.br.'


* Did you recover the Wuhan-1 genome?
* How similar is it to the reference?
* What do all of the fields mean (use <code>blastn -help</code> to see the full list of blast options)?



