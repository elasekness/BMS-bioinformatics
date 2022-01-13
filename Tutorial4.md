## Tutorial 4

**Objective:** To Perform an RNA-Seq analysis and become familiar with the basics of R.

RNA-Seq analyses are employed to study the effects of a treatment (such as antibiotic exposure) on the expression levels of
an organism's or cell's transcriptome.  When a reference sequence is available, we can map our reads to it and obtain count data.
These read counts serve as a proxy of expression for genomic elements of interest (i.e. genes). When no reference assembly is available,
a de novo assembly of the transcriptome is performed, to which reads are mapped back. Thus, one of the first steps is deciding which
type of assembly you will be performing.  Another essential consideration is the experimental design.  You will need a control against which
your treatment data can be compared and you will need biological replicates for both treatments and controls. Ideally, only the variable being
tested should vary between the biological replicates for your treatment and the control conditions should be kept as similar as possible. The number of
replicates is also important and most agree that at least three should be used.  A good reference regarding the number of biological replicates
to consider is [Schurch et al 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4878611/). Finally, you will need
to decide which tool to employ for estimating differential expression.  As with assemblers, there are many tools available for performing differential
expression analyses.

<br>

## Download and install R on your computer

If you don't have R on your computer, please install it: [https://www.r-project.org/](https://www.r-project.org/)

R is freely available software for statistical analyses and figure generation.  It is one of the main tools employed by bioinformaticians.
R packages are bundles of code that perform myriad operations and analyses and can be installed from a CRAN mirror.  Bioconductor is essentially
a project for developing and maintaining code, written in R, for analysis of biological data.
We will install [Bioconductor](https://www.bioconductor.org/)
as well as one of its packages: [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).

Open an R window on your computer by double-clicking on the application and follow the instructions from the links provided above to install Bioconductor and DESeq2.

<br>


## Download the fastq files for our RNA-Seq analysis from the SRA

We will be analyzing a dataset of <i>Staphylococcus aureus</i> N315 cultures exposed to the antibiotic daptomycin at 4 ug/ml.
Daptomycin is a cyclic lipopeptide which binds to the cytoplasmic membrane of <i>S. aureus</i> in a calcium-dependent manner,
ultimately leading to membrane depolarization and cell death. Daptomycin non-susceptibility is most frequently associated
with mutations in the multi-peptide resistance factor gene (mprF). The mprF gene encodes for lysyl-phosphatidylglycerol (LPG) synthetase,
responsible for lysinilating phosphatidylglycerol (PG) and translocating it to the outer membrane.
Mutations in mprF are thought to mediate daptomycin resistance by causing an increase in the LPG content of the cell membrane,
which increases the net positive charge of the membrane and repels daptomycin.
Mutations in the cardiolipin synthase gene (cls2) have also been observed in DAP-NS strains
and might further alter membrane composition by increasing the ratio of cardiolipin to PG, thereby leaving fewer target sites for daptomycin.
Other phenotypic changes that often appear in DAP-NS strains include increased cell wall thickening and changes in membrane fluidity
that might prevent the antibiotic from accessing its target site. Thus we might also expect changes in gene expression for genes involved in 
cell wall and/or the cell membrane homeostasis.

There are two biological replicates for the treatment and the control, which you can access via the bioproject [PRJNA669520](https://www.ncbi.nlm.nih.gov/bioproject/?term=669520)
by clicking on the SRA experiments link.

Select the four experiments of interest
* 4 ug/ml daptomycin treated samle, grown in minimal media (SSM9PR)
* 4 ug/ml daptomycin treated samle, grown in minimal media (SSM9PR)
* WT control treated with DMSO and grown in minimal media (SSM9PR)
* WT control treated with DMSO and grown in minimal media (SSM9PR)

and send the results to the 'Run Selector' to view the four SRR accession numbers in one location.

| Run | Library name |
| --- | ------------ |
| SRR12830230 | Daptomycin-2 |
| SRR12830233 | UT-control-2 |
| SRR12830234 | Daptomycin-1 |
| SRR12830237 | UT-control-1 |

Make a new directory for this analysis, `cd` into it and download the fastq files for each run to your VM with `fasterq-dump`.

	mkdir rnaseq
	cd rnaseq
	fasterq-dump
	fasterq-dump SRR12830230
	fasterq-dump SRR12830233
	fasterq-dump SRR12830234
	fasterq-dump SRR12830237

Rename the PE fastq files something meaningful.

	mv SRR12830230_1.fastq dap2_1.fastq
	mv SRR12830230_2.fastq dap2_2.fastq
	mv SRR12830233_1.fastq cont2_1.fastq
	mv SRR12830233_2.fastq cont2_2.fastq
	mv SRR12830234_1.fastq dap1_1.fastq
	mv SRR12830234_2.fastq dap1_2.fastq
	mv SRR12830237_1.fastq cont1_1.fastq
	mv SRR12830237_2.fastq cont1_2.fastq

> **Note:** The renamed fastqs are also available in our GCP bucket: gs//wc-bms-bi-training-bucket/rnaseq/fastq

<br>

## Clean your reads with TrimGalore

These reads were generated on an Illumina NextSeq in a 2x50 bp format.
You can process your PE fastqs one-by-one or you could execute a BASH for-loop to do the job for you.

The long way:

	trim_galore -q 30 --paired cont1_1.fastq cont1_2.fastq
	trim_galore -q 30 --paired cont2_1.fastq cont2_2.fastq
	trim_galore -q 30 --paired dap1_1.fastq dap1_2.fastq
	trim_galore -q 30 --paired dap2_1.fastq dap2_2.fastq

The short way:

	ls *_1.fastq | cut -d "_" -f 1 > seqlist
	for filn in `cat seqlist`; do trim_galore -q 30 --paired $filn"_1.fastq" $filn"_2.fastq"; done

> Here, we are listing all R1 fastq files and cutting them on the underscore delimiter to take the first field, which is the
> base name for each PE fastq file. We then save that output to a file called 'seqlist.'  You can `cat` the seqlist file to see the results.
> The next command is the for loop.  Instead of looping through files one-by-one, we are looping through the 'seqlist' file line-by-line
> to obtain the basename of each PE run. The back ticks represent a subprocess. The output of the subprocess command `cat` is being passed
> to the for loop.  Thus, each basename in our 'seqlist' file becomes a variable.  We then use this basename to specify the R1 and R2 fastq files
> by filling in the reminder of the unique part of each PE file's name.  For example, `$filn` will get interpreted as 'cont1' and the file endings
> "\_1.fastq" and "\_2.fastq" will be interpreted literally because of the quotation marks so that we get the full file names, 
> 'cont1_1.fastq' and 'cont1_2.fastq'.

<br>

## Download a reference genome for _S. aureus_ N315

Return to NCBI and search 'assembly' for 'Staphylococcus aureus N315.'

> **Note** We previously searched genomes for SARS-CoV-2 and browsed a list of genomes to access the Wuhan-1 reference RefSeq ftp site.  Here we are taking an 
> alternative route.

Click on the link to the RefSeq FTP page, which is to the right on the for the _S. aureus_ N315 assembly page.

Copy the genome assembly and annotation files to your VM. We'll also want the coding sequences for an alternative approach to read quantification.

	curl -O "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/645/GCF_000009645.1_ASM964v1/GCF_000009645.1_ASM964v1_genomic.fna.gz"
	curl -O "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/645/GCF_000009645.1_ASM964v1/GCF_000009645.1_ASM964v1_genomic.gff.gz"
	curl -O "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/645/GCF_000009645.1_ASM964v1/GCF_000009645.1_ASM964v1_cds_from_genomic.fna.gz"

<br>

## Map your reads to the reference genome for _S. aureus_ N315

	gzip -d *gz
	cp GCF_000009645.1_ASM964v1_genomic.fna sa.fna
	bwa index sa.fna 
	for filn in `cat seqlist`; do bwa mem sa.fna $filn"_1_val_1.fq" $filn"_2_val_2.fq" | samtools sort | samtools view -F 4 -o $filn".sorted.bam"; done
	
> We are decompressing our assebmly and annotation file, copying the assembly file to a shorter name, indexing it for `bwa`, running `bwa` to map our reads, 
> then using `samtools` to sort our alignment and convert it to bam format, while only including reads that aligned to our reference genome.
> **Note:** BAM alignments are also located in our GCP bucket: wc-bms-bi-training-bucket/rnaseq/bams

<br>

## Count the number of reads that overlap coding sequences

We can use the function `multiBamCov` from the package [bedtools](https://github.com/arq5x/bedtools2) to count the number of reads that overlap genic regions. The `multiBamCov` function requires the bam file as well as the coordinates of your genes in a bed or gff format.

First, we'll need to index our sorted bam files. You can use a for loop to automate the process.

	for filn in *sorted.bam; do samtools index $filn; done

In this experiment, we are only interested in the expression of genes but our gff file contains annotation for genes and their coding sequences, as well as other genomic elements (tRNAs, rRNAs, pseudogenes). Therefore it would be easier to convert our gff file to a bed file that only contains the genomic intervals of interest.

	grep -P "\tCDS\t" GCF_000009645.1_ASM964v1_genomic.gff | cut -f 1,4,5 | grep -v "NC_003140" > sa.bed
	
> The `P` argument indicates the use of regular expressions in your `grep` command.  Here we are searching for all lines that have 'CDS' pre- and pro-ceeded by 
> a tab. We are then cutting the first,fourth, and fifth fields of the output, which correspond to the genome name and the coding sequence coordinates.  
> The _S. aureus_ genome also contains a plasmid.  Our final `grep` command specifies to exclude any lines that contain the plasmid accession number.

Perform the read count, finally.

	multiBamCov -bams cont1.sorted.bam cont2.sorted.bam dap1.sorted.bam dap2.sorted.bam -bed sa.bed > sa-bwa.counts.txt

> The output is saved as a tab-delimited file, where each line contains information for our genomic intervals of interest (CDSs).
> The first field is the genome accession number, followed by the CDs coordinates, and the read counts for the bam files in the order
> you specified in your `multiBamCov` command.
> With a few adjustments, the output of this command can be used as the count matrix for our DESeq2 analysis.

Convert the genome accession numbers on each line back to the coding sequence name (names must be unique for downstream analyses) and
elminate the CDS coordinates.

	grep -P "\tCDS\t" GCF_000009645.1_ASM964v1_genomic.gff | grep -v "NC_003140" | cut -f 9 | cut -d "=" -f 3 | sed "s/gene-SA_//" | sed "s/;.*//" > 		sa_tags.txt

> This command is parsing our GFF file to return only the gene name or locus tag associated with our coding sequences. 
> However, the locus tag 'RS04035' is actually listed twice in our file.  We can edit one of the duplicate names in nano so that it is unique.
> `Ctrl-w` in nano allows you to perform a search.

	paste sa_tags.txt sa-bwa.counts.txt | cut -f 1,5-8 > sa-bwa.countsR.txt 

> This command is attaching our locus tag names to the read count data and eliminating the fields we don't need for our DESeq2 analysis.
> `paste` pastes two files side-by-side with a tab in between them.
> Add headers to your 'sa-bwa.countsR.txt' file and it's ready for import into R!
> **Note:** 'sa-bwa.countsR.txt' is also in our GCP bucket: gs://wc-bms-bi-training-bucket/rnaseq/readcounts


Download the read count file to your computer if you have R installed.  You can also use the command-line version of R on our VMs.

<br>

## Follow the demo in class for a brief introduction to R and DESeq2

We will now perform our differential expression analyses in R using DESeq2. There are several ways to import data into DESeq2 but we will follow the 
instructions for importing a count matrix. After we have imported are data into R, the basic commands will be:

	R
	library(DESeq2)
	countsTable = read.table("sa-bwa.countsR.txt", header=T, sep="\t", row.names=1)
	ColData = read.table("ColData.txt", header=T, sep="\t", row.names=1)
	dds <- DESeqDataSetFromMatrix(countData = countsTable, colData = ColData, design = ~ Condition)
	dds <- dds[rowSums(counts(dds)) >100, ]
	dds <- DESeq(dds)
	res <- results(dds, contrast=c("Condition", "treated", "control"))
	write.table(dds, file='sa-bwa_deseq.txt', sep='\t')

> Here we are importing both our read count table and a table called ColData, which we must create.  This contains information on which columns in our read
> count table are controls and which are treatments. Both of these are used to create the DESeq2 dataset 'dds.' Note that un-normalized read counts 
> should be used as DESeq2 will correct for differences in library size. We also need to specify a design, or which variables to fit in our model and test for
> differential expression (DE). Here we are testing the 'Condition' (treated vs. control) on the effects of gene expression.
> The DESeq command will be discussed more in class.
> If you're working on your computer, you can make a 'ColData.txt' file in Excel.  If you're working on our VMs, you can make the same table in `nano`.

<br>

## Quantify transcripts with Salmon

[Salmon](https://combine-lab.github.io/salmon/) is a 'wicked fast' program that allows the direct quantification of reads against a transcriptome (no need for SAM or BAM files).  Because we don't have a transcriptome (we did not perform a de-novo assembly of our reads), we'll use the coding sequences that we downloaded as a substitute.

Copy the cds file and name it something simple.

	cp GCF_000009645.1_ASM964v1_cds_from_genomic.fna sa_cds.fna

If you look at the definition lines for the coding sequences, you'll see that the names are quite long and start with 'lcl|.' Let's remove that from all the definition lines to make things a bit cleaner.

	sed "s/lcl|//" sa_cds.fna | grep ">" | head
	sed -i "s/lcl|//" sa_cds.fna
	
> First check to see that the substitution command worked as you intended.
> Then make the replacement. The `-i` indicates that sed should make the substitutions in place.

Index the coding sequence file.

	salmon index -t sa_cds.fna -i salmon_index
	
> If you `ls` your directory, you'll see that salmon has put the index files to your 'transcriptome' in a sub-directory called 'salmon_index'

Quantify transcript expression.

	for filn in `cat seqlist`; do salmon quant -i salmon_index -l A -1 $filn"_1_val_1.fq" -2 $filn"_2_val_2.fq" --validateMappings -o "salmon_quant/"$filn; done
	
> The `quant` command allows direct quantification of reads agains the 'transcriptome' index and will output the results into a directory called 'salmon_quant.'
> The TSV files for each library will be in a subdirectory named by the basename of the file.  The `-l` option specifies the automatic detection of the 
> library type.  The `--validateMappings` option is the recommended default.  It essentially checks that the mappings are plausible enough to be quantified.

* Was this considerably faster and simpler than aligning with BWA?
* Notice how we don't need additional steps to generate a BAM file and then extrat read counts from it.
* Look at the TSV files for your libraries (basename.sf).  You are given both TPM (Transcripts Per Million) and raw read count quantifications.
* The multiBamCov command allows you to put readcounts from multiple bam files into the same file but Salmon outputs these to separate files.

## Perform a DESeq2 analysis with your Salmon-quantified reads.

* Are the number of DEGs the same between the BWA method and Salmon?
* Did they find similar read counts for your coding sequences?
* Are the DE genesets the same?
* Which method do you prefer?
	

