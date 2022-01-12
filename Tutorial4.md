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
> by filling in the reminder of the unique part of each PE file's name.  For example, `$filn` will get interpreted as `cont1` and the file endings
> "\_1.fastq" and "\_2.fastq" will be interpreted literally because of the quotation marks so that we get the full file names, 
> 'cont1_1.fastq' and 'cont1_2.fastq'.


