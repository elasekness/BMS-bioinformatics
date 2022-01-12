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

**Note:** Now that our VMs have additional RAM, we should be able to use all the reads. Below, we are running the SPAdes installed on our VM and foregoing the read error correction step to save time.

	spades --only-assembler -1 SRR10971381-trimmed_1.fq -2 SRR10971381-trimmed_2.fq -o wuhan_assembly_spades

> SPAdes provides several assembly options (including the metagenome option) and is ideally suited for smaller genomes.
> SPAdes will perform error-correction prior to assembling or you can run without the error-correction step as we have done in the second example.
>
> **Note:** The ampersand specifies performing a job in the background. This returns your command-line prompt so that you can continue working.
> Given the size of the library (16 GB!), assembly could take quite some time if we used all reads. You can see that your job is stil running by typing `top`,
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

<br>

## Assemble the same reads with Shovill

	docker run --rm -v $(pwd):/data -w /data staphb/shovill shovill --R1 SRR10971381_1.fastq --r2 SRR10971381_2.fastq --trim --outdir wuhan_assembly_shovill

> As mentioned previously, Shovill is a faster implementation of SPAdes. Notice how much faster Shovill ran even with all of our data!
> We are using the `--trim` option to remove adapters but we could have started with out trimmed reads.
> The one caveat is that Shovill doesn't have a metagenome mode.

Briefly explore the Shovill assembly and its summary output to see how the program compares to SPAdes.

The final contigs are in a file called 'contigs.fa' and are conveniently organized by descending length.
Also notice that the naming conventions for the contigs are much more manageable (although the original name assigned by SPAdes is also included).

* Did Shovill recover the Wuhan genome?
* What assembler would you choose based on this exercise?



## Perform a metagenome assembly with MEGAHIT


	megahit -1 SRR10971381-trimmed_1.fastq -2 SRR10971381-trimmed_2.fastq -o wuhan_assembly_megahit

> Notice how much faster MEGAHIT completed in comparison to SPAdes.
> Your contigs are located in the file called 'final.contigs.fa' in your output directory 'wuhan_assembly_megahit'.

* Generate the same summary statistics for your MEGAHIT assembly as you did for your SPAdes assembly.
* Did all three assemblies recover the Wuhan-1 genome?
* If so, are they similar lengths and coverage?
* Other tools to try are: [Pangolin](https://pangolin.cog-uk.io/) and [Nextclade](https://clades.nextstrain.org/).

<br>


## Compare the performance of SPAdes, Shovill, and Megahit with [QUAST](https://github.com/ablab/quast)


Quast is a genome (and metagenome) evaluation tool that can compare the quality of assemblies produced by different programs.
It will output the number of large contigs for each assembly, the length of the longest one, the N50, and the number of genes predicted.
When a reference assembly is provided, Quast will output additional statistics, such as the number of misassemblies, the percentage of the reference genome
recovered, etc.

Copy your contig files from your three assembly directories to the directory with your reference assembly and rename them something meaningful


	cp ~/wuhan_fastqs/wuhan_assembly_spades/scaffolds.fasta ~/reference_assembly/spades.fa
	cp ~/wuhan_fastqs/wuhan_assembly_shovill/contigs.fa ~/reference_assembly/shovill.fa
	cp ~/wuhan_fastqs/wuhan_assembly_megahit ~/reference_assembly/megahit.fa

> This command assumes your assembly output directories are in a directory called 'wuhan_fastqs' and the reference assembly and annotation file for Wuhan-1
> from NCBI is in the directory 'reference_assembly.' The `~` symbolizes your home directory.

<br>

Run the Dockerized version of Quast to comare the three assemblies.

	cd wuhan_fastqs
	docker run --rm -v $(pwd):/data -w /data staphb/quast quast.py -r wuhan.fna -g GCF_009858895.2_ASM985889v3_genomic.gff spades.fa shovill.fa megahit.fa

> Your results will be located in a directory called 'quast_results' in a subdirectory called 'latest.' 
> Notice the summary report is given to you in multiple formats.

<br>

Download the html version of the results ('report.html') to your computer.

The html report provides a nice visual of the results in table and figure formats. Because
we provided a reference genome and a GFF annotation file, Quast compares the performance of
all three assemblers in recovering a complete and accurate (misassemblies, INDELs, etc) Wuhan-1 genome
as well as the number of complete genomic features.

* Which assembler performed best in assembling the Wuhan-1 genome?
* What characteristics would you want to see in a metagenome assembly?
* Given those characteristics, which created a better metagenome assembly?

<br>

## Annotate your consensus SARS-CoV-2 genome with [Prokka](https://github.com/tseemann/prokka)

We will use a program called 'Prokka' to annotate our consensus genomes that we generated in Tutorial 2. Prokka takes advantage of several other tools to make high-quality predictions for coding sequences, tRNAs, rRNAs, and CRISPRs. Prokka uses BLAST against well characterized protein databases as well as HMM (Hiden Markov Models) scans. Hidden Markov Models are probabilistic models about the identity of an amino acid/nucleotide at a certain position in a protein/gene generated from multi-sequence alignments.
HMM based predictions can be particularly powerful in identifying distantly related homologs that a Blast search wouldn’t detect.
<br>

Run the dockerized version of Prokka in interactive mode to see the help menu.

	docker run -it staphb/prokka

> As you can see, there are many options associated with Prokka, including ways to customize the annotation output.
> Some useful arguments include specifying the name of the locus tag, adding the genus and species names to the annotation files,
> naming the output directory, and specifying the annotation mode - in this case it will be 'Viruses.'

<br>

Before running Prokka, change the ownership of your consensus genome fasta file and shorten the sequence name.

	sudo chown root:your_username IDRnumber.fa

> `chown` = change ownership. Since we generated our consensus genomes with a dockerized version of iVAR without changing the user, the fasta file belongs to 
> 'root.' This means that we cannot modify the file.  We need to shorten the name of the definition line (or sequence name) because Prokkka will complain about
>  anything longer than 20 characters.  Thus we must change the ownership of the file from root to us.  After you have changed the ownership of this file, you can 
>  shorten the sequence name in Nano.

Instead of changing the ownership of the fasta file, you could also change the name of the sequence with `sed` and write the output to a new file. You can first use `grep` to obtain the name of the sequence.

	sed 's/Consensus_\(.*)_S.*/\1/' IDRnumber.fa > IDRnumber_modified.fa

> We'll talk about the syntax of this command in class. The `\1` and escaped parentheses have special meaning.

## Run Prokka with the 'Viruses' option for annotation mode

	docker run --rm -v $(pwd):/data -w /data staphb/prokka prokka --kingdom 'Viruses' --species 'SARS-CoV-2' --prefix 'IDRnumber' --outdir IDRnumber_annotation IDRnumber.fa

* Look at the output generated.
* How many proteins were annotated?
* Do those seem correct? How could you test this?
* Do the coding sequences have the same coordinates as the reference genome? Hint: compare GFF files.

## If there's time, annotate with [VADR](https://github.com/ncbi/vadr/wiki/Coronavirus-annotation#howto)

As you can see, there are many ways to get the same answers!








