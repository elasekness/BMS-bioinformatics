## Tutorial 1

**Objective:** Successfully log on to Google Cloud Platform (GCP) Console and become familiar with working in a command-line (Linux) environment.

<br>

## Opening a terminal window

We are working on remote servers or virtual machines (VMs) hosted by Google.
To open a terminal window where we will work, we need to log on to the GCP console with our Wadsworth credentials,
start the VM associated with our Project and `ssh` to those servers.
`ssh` is a secure shell protocol for safely connecting to remote services.
To `ssh` you need an account on the server with a username and password as well as the IP address of the server.

- Navigate to the [GCP console](https://console.cloud.google.com) and log on with your Wadsworth credentials
- Once connected to the console, click the Compute Engine link and start your VM by expanding the three dot icon to the right
- Click the SSH button once the VM has started. This will open a terminal window on your VM.

To perform any operations in a Linux environment, we need to tell the computer what to do by typing specific commands into our terminal window.
In general, the syntax will be: `Command File`,
where `Command` is the function you want to perform on the` File` or `Directory`
that you specify. Sometimes specifying a command is all we need.

Examine the contents of your directory with:
 
`ls`

 > **`ls`** = list

<br>

## Directory structure and navigation

Directories are folders with specific locations on the server.
To navigate from one directory to another, we must specify the location of the directory or its path.
Directories have a forward slash `/` after their names.  Thus to get to a subdirectory within a directory
you would specify the path by stringing together the directory names, separated by `/`. We'll practice below

First, determine the name and location of your home directory.

`echo $HOME`

> **`echo`** = repeat, **`$HOME`** = variable name for your home directory and its location on the server.
> Try the `echo` command with other variable names, such as **`$SHELL`** or **`$PATH`**.

<br>

You can also print your working directory (where you are).
 
`pwd`

> **`pwd`** = print working directory

<br>

Make a directory within your home directory called `genomes`.

`mkdir genomes`

> **`mkdir`** = make directory

<br>

Change (move) to a different directory.

`cd genomes`

> **`cd`** = change directory. Use **`../`** or **`..`** to move up one directory (back to your home directory). What does **`cd`** alone do?

<br>


## Permissions

Directories and files have specific permissions associated with them or things that you are allowed to do to them.
There are three permissions: 1) the ability to read (r) 2) the ability to write (w) and 3) the ability to execute (a script or program; x).
There three sets of permissions representing 1) the User (you), 2) the group (multiple users),
and 3) Others (Everyone else in the world).
<br>

Examine the permissions of your 'genomes' directory.

`ls -l genomes`

> **`ls -l`** = long list, providing you information on when the **`genomes`** directory was created and its associated permissions.

<br>

Change the permissions associated with your `genomes` directory with the 'chmod' command (change mode).

`chmod 775 genomes`

> Permissions are represented by three-digit (for user, group, and other) octal numbers.
> Here we are allowing the user and group universal permissions (7 = read, write, and execute) and all others
> the ability to read and write only (5).

> For more information on permissions, see: [Linux permissions](https://www.guru99.com/file-permissions.html#linux_file_ownership)

<br>


## Manipulating files (making, (re)moving, editing, and viewing)

There are many ways to make and view files in a Linux OS (operating system).
We can redirect output from a command that prints its output to your screen (STDOUT) to a file instead,
we can generate files on the fly by opening them in a text editor, or we can copy an existing file.
Similarly, we can view and edit files in a text editor or we can print their contents to the screen with various command-line options.
As a general rule, it's always good to examine some of the contents of your file to ensure you've generated the results you want in the
format you want it. Or that you are using the correct file and format for downstream applications.

Redirect STDOUT to a file.

`ls /usr/bin > programs.txt`

> The path **`/usr/bin`** specifies the location where various Bash commands are found. When you type a command, **`/usr/bin`** is one of the locations
> your computer searches to find and execute the command. Was **`/usr/bin`** part of your **`$PATH`**?
> Here we are redirecting the STDOUT from the **`ls`** command to a file named **`programs.txt`**. The **`>`** sign is responsible for the redirection.

<br>

Scroll through the contents of your file.

`more programs.txt`

> Scroll through line-by-line with the enter key.  Scroll through page-by-page with the space key.
> Do you notice that the file contains some of the commands you have just used?
> Exit with **`control-c`**.

<br>

Display the first ten lines of your file.

`head programs.txt`

> **`head`** displays the first ten lines by default but you can specify the number of lines with a flag.
> For example, **`head -200 programs.txt`** will display the first two hundred lines of your file. In general, most
> commands have additional arguments (or flags) associated with them. You can see the different usage statements
> by typing the command with a **`-help`** or **`--help`** option.

<br>

Display the last ten lines of your file.

`tail programs.txt`

<br>

Print the entire contents of your file to your screen.


`cat programs.txt`

> **`cat`** = concatenate.  The **`cat`** command can also join the contents of multiple files together.

<br>

Make and view the contents of a file with a text editor.


`nano new_file.txt`

> Nano, emacs, vim, and vi are all text editors.
> You can make an empty file on the fly as we did here.  This will open a blank text editor screen.
> Type some content and save it with **`control-o`**. To exit the text editor, use **`control-x`**.
> More information on Nano commands can be found here: 
> [Nano](https://www.howtogeek.com/howto/42980/the-beginners-guide-to-nano-the-linux-command-line-text-editor/)

<br>

Rename a file.

`mv programs.txt installed_programs.txt`

> **`mv`** = move. Renaming files with `mv` will overwrite any existing file.  You can also mv a file to a different directory.  
> Try it: **`mv installed_programs.txt genomes`**
> Can you move the file back to your home directory?

<br>

Copy a file.

`cp installed_programs.txt duplicate.txt`

> **`cp`** = copy. Can you copy a file to a different directory in one command?

<br>

Remove a file.

`rm duplicate.txt`

> **`rm`** = remove.  Remember, a removed file cannot be restored.  
> Can you remove a file from a different directory without having to change directories? 
> How would you remove a directory?

<br>

## Notes on working in a Linux environment


* Avoid creating file names and directories with spaces. Use underscores, dashes, or periods instead to separate multi-part names.
* Spaces need to be escaped in Linux (more on that later).  For example, if you tried to make a directory called “my directory”,
* mkdir would make two directories, “my” and “directory.”
* Use autocomplete for speedier typing and to avoid typos.  Autocomplete will fill in the unique part of a command or file.
For example, if I had only one file in my directory that began with a “b,” I could type b and then press the tab key to autocomplete the name of the directory.
* Everything in Linux is case sensitive
* Can’t find a command?
Try `which` to see if the command is in your path – whether the command is in a location that the computer searches for executing commands.
* Hit the up-arrow key to recall a command you entered previously.

<br>

## Databases and obtaining sequences


There are several sequence databases – NCBI, JGI, EMBL - that you might encounter.
You might want to explore each of these to familiarize yourself with the resources they offer.
We will focus on NCBI.  Our goal is to download a reference genome (Wuhan-1) for SARS-CoV-2,
which we will use later to perform a reference-based assembly and an outbreak analysis.
We will also want to download the accompanying annotation file (gff file), which provides a description of the genes,
the function of the coding sequences, and the nucleotide positions of the genes in the genome and the translated coding sequences (faa file).


Navigate to NCBI’s homepage: [https://www.ncbi.nlm.nih.gov/](https://www.ncbi.nlm.nih.gov/)


> Notice that there are options to submit sequences, download sequences, and even analyze data.
> PubMed allows literature searches and BLAST is an alignment tool to look for sequences that are similar (a proxy for homology) to your queries.
> Also notice that NCBI has provided a quick link to SARS-CoV-2 data. You could obtain a nucleotide record by clinking on this link but we'll follow
> the more traditional route for now.


Choose “Genome” under the top left pull-down menu (set to “All Databases” by default), type SARS-CoV-2 into the search area, and hit enter.
This brings us to a page containing information on the reference genome (Wuhan-1).  We could also 'Browse the list' of other available genomes.
We can click on the 'RefSeq' link and use the 'Send to' menu to save the 'Complete Record' for this genome to a file in 'fasta' format. By default,
the file is named 'sequence.txt.'  From there we can upload the file to our VM. Instead we will:


'Browse the list' of available genomes to access the RefSeq or GenBank FTP site for our reference genome (which is conveniently the first genome listed).
 If it wasn't the first genome listed, we could apply search filters to narrow the list.


Click on the RefSeq ('R') FTP link for MN908497.3 (Wuhan-1). Open the link as a "Guest." This takes you to a directory with several files associated with the Wuhan-1 genome. For now let's obtain the genome assembly (.fna), the accompanying annotation (.gff), and the protein coding sequences (.faa).
Instead of a two-step process of downloading files to our computers and then uploading them to the cloud, we will use the **`curl`** command (copy url) to copy the files directly to
our VMs.


Right click on the fna file and 'Get info.' The information page will list the file location on NCBI's server.


Copy the server information and type the following command in your VM terminal:

`curl -O "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz"`

> The **`-O`** option saves the copied contents to a file named as it was on the FTP site.
> Repeat this for the gff and faa files

<br>

## Obtaining reads from the SRA


Now let's download the raw reads for the Wuhan-1 reference genome from the SRA (sequence read archive) database.  We'll use these later to perform a de-novo assembly. Typically, any published NGS data must also be submitted to the SRA. Each sample/specimen sequenced will have a BioSample accession number. Biosample information provides associated metadata. The SRA and Biosample for each submission are further housed under a BioProject, which can contain multiple submissions from the same study or experiment


Return to the [NCBI homepage](https://www.ncbi.nlm.nih.gov/) and search genomes for Wuhan-1 again. This brings you to its GenBank page.
Normally, we would click the link to the associated BioProject but this takes us to a large umbrella project that would be difficult to search.
Instead, we could click the link to the publication on PubMed and search the journal article for the BioProject.  I have done this for you.


Go to the Bioproject for Wuhan 1: [PRJNA603194](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA603194)
Notice the BioProject page gives you information on the purpose of the study, as well as several links, including one to the SRA experiment.


Click on the link to the [SRA experiment](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=603194) for Wuhan-1
Notice this is a metatranscriptome, meaning this may not be a pure isolate of a single SARS-CoV-2 virus (although human RNA/DNA should have been removed prior to
SRA submission).  The sequencing was performed on an Illumina MiniSeq and the reads are paired (i.e. there are Forward and Reverse reads for the same amplicon).
We could download the two fastq files using the NCBI link provided or we could use faster tools provided by NCBI.


Use **`prefetch`** and **`fasterq-dump`** tools from the SRA toolkit to download the Wuhan-1 fastq files.

`prefetch SRR10971381`
`fasterq-dump SRR10971381`

> **`prefetch`** will download the SRA data in binary format and **`fasterq-dump`** will perform the fastq conversion.
> Notice that the conversion tool automatically saves forward and reverse reads to separate files.
> Although we could use fasterq by itself, NCBI claims prefetch in combination with fasterq is faster.
> More information on the SRA-toolkit and other frequently used tools can be found here: [SRA toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc)
>
> **Note:** If your are having trouble obtaining the fna, faa, and gff files from NCBI, you can also copy them to your account on the VM from our GCP bucket
> using the gsutil commands provided by GCP. Fastq files also available in the bucket in the same folder: 'reference_assembly.' However, you will need to
> authenticate your account first in order to use the gsutil commands.
>
> `gsutil -m cp gs://wc-bms-bi-training-bucket/reference_assembly/GCF* .`
>
> We'll discuss this command more in class but essentially, we are copying all files that start with GCF from our GCP bucket to our current directory `.` .


