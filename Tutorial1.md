## Tutorial 1

Objective: Successfully log on to Google Cloud Platform (GCP) Console and become familiar with working in a command-line (Linux) environment.
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

> **`echo` ** = repeat, **`$HOME`** = variable name for your home directory and its location on the server.
> Try the `echo` command with other variable names, such as **`$SHELL`** or **`$PATH`**.

<br>

You can also print your working directory (where you are).
 
 `pwd`

 `pwd` = print working directory

<br>

- Make a directory within your home directory called `genomes`.

 `mkdir genomes`

 `mkdir` = make directory

<br>

- Change (move) to a different directory.

 `cd genomes`

 `cd` = change directory. Use `../` or `..` to move up one directory (back to your home directory). What does `cd` alone do?

<br>


## Permissions

Directories and files have specific permissions associated with them or things that you are allowed to do to them.
There are three permissions: 1) the ability to read (r) 2) the ability to write (w) and 3) the ability to execute (a script or program; x).
There three sets of permissions representing 1) the User (you), 2) the group (multiple users),
and 3) Others (Everyone else in the world).
<br>

- Examine the permissions of your 'genomes' directory.

 `ls -l genomes`

 `ls -l` = long list, providing you information on when the `genomes` directory was created and its associated permissions.

<br>

- Change the permissions associated with your `genomes` directory with the 'chmod' command (change mode).

 `chmod 775< genomes`

Permissions are represented by three-digit (for user, group, and other) octal numbers.
Here we are allowing the user and group universal permissions (7 = read, write, and execute) and all others
the ability to read and write only (5).

For more information on permissions, see: [Linux permissions](https://www.guru99.com/file-permissions.html#linux_file_ownership)

<br>


## Manipulating files (making, (re)moving, editing, and viewing)

There are many ways to make and view files in a Linux OS (operating system).
We can redirect output from a command that prints its output to your screen (STDOUT) to a file instead,
we can generate files on the fly by opening them in a text editor, or we can copy an existing file.
Similarly, we can view and edit files in a text editor or we can print their contents to the screen with various command-line options.
As a general rule, it's always good to examine some of the contents of your file to ensure you've generated the results you want in the
format you want it. Or that you are using the correct file and format for downstream applications.

- Redirect STDOUT to a file.

 `ls /usr/bin > programs.txt`

The path `/usr/bin` specifies the location where various Bash commands are found. When you type a command, `/usr/bin` is one of the locations
your computer searches to find and execute the command. Was `/usr/bin` part of your `$PATH`?
Here we are redirecting the STDOUT from the `ls` command to a file named `programs.txt`. The `>` sign is responsible for the redirection.

<br>

- Scroll through the contents of your file.

 `more programs.txt`

Scroll through line-by-line with the enter key.  Scroll through page-by-page with the space key.
Do you notice that the file contains some of the commands you have just used?
Exit with `control-c`.

<br>

- Display the first ten lines of your file.

 `head programs.txt`

`head` displays the first ten lines by default but you can specify the number of lines with a flag.
For example, `head -200 programs.txt` will display the first two hundred lines of your file. In general, most
commands have additional arguments (or flags) associated with them. You can see the different usage statements
by typing the command with a `-help` or `--help` option.

<br>

- Display the last ten lines of your file.


 `tail programs.txt`


<br>

- Print the entire contents of your file to your screen.


`cat programs.txt`

 `cat` = concatenate.  The `cat` command can also join the contents of multiple files together.

<br>

- Make and view the contents of a file with a text editor.


 `nano new_file.txt`

Nano, emacs, vim, and vi are all text editors.
You can make an empty file on the fly as we did here.  This will open a blank text editor screen.
Type some content and save it with `control-o`. To exit the text editor, use `control-x`.
 More information on Nano commands can be found here: [Nano](https://www.howtogeek.com/howto/42980/the-beginners-guide-to-nano-the-linux-command-line-text-editor/)

<br>

- Rename a file.


 `mv programs.txt installed_programs.txt`

`mv` = move. Renaming files with `mv` will overwrite any existing file.  You can also mv a file
to a different directory.  
Try it: `mv installed_programs.txt genomes`
Can you move the file back to your home directory?

<br>

- Copy a file.


 `cp installed_programs.txt duplicate.txt`

<code>cp</code> = copy. Can you copy a file to a different directory in one command?

<br>

- Remove a file.


 `rm duplicate.txt`

<code>rm</code> = remove.  Remember, a removed file cannot be restored.  Can you remove
a file from a different directory without having to change directories? How would you remove a directory?

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


