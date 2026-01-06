### Sequence Bioinformatics - Group Project - WS 2025
# Server Instructions
 
## access the server:

`ssh username@sshgw.cs.uni-tuebingen.de`

Once you're in the gateway, you can type (for one of the two servers):

`ssh nog`

or

`ssh collins`

## Please work in this folder:

`/teachstor/share/groupprojectWS25/groupB/`

## BEFORE starting to work

The following commands need to be executed, before starting any pipeline on the server.
You need the `bash` command to start conda service. I put them here to copy/paste them.

```
ssh [USERNAME]@sshgw.cs.uni-tuebingen.de
ssh nog
cd ../../teachstor/share/groupprojectWS25/groupB/
bash
conda activate metagenomics
```

## AFTER work: change permissions of files
after uploading or creating files (through running the pipeline etc.) do:

```
cd teachstor/share/groupprojectWS25/groupB
chmod -R 770 ./*
```

## miniconda

You do not have sudo rights, so you cannot install anything using `sudo apt install`. You always need to either download some installer/binary or install packages using conda.

To install miniconda, you can follow these instructions:
```
https://www.anaconda.com/docs/getting-started/miniconda/install#linux-terminal-installer
```

(open "macOS/Linux installation" dropdown and go to "Linux terminal installer").

<b>
conda config in .bashrc might have to be copied to everyones home directory.
 If init of conda does not happen automatically upon entering the server type:
</b>

`bash`
## linux version check
`uname -m`


## install tools
If there are any tools you need to install that you cannot install using conda
(you can always google "conda install tool-name" to find that out),
the tools usually provide instructions on how to install them via the command line.
You can add them to your PATH variable if you want to.
For this, add something like this to your `.bash_profile`
(if you do not have that file, you can create it in your home directory):
```
export TOOL_HOME=/teachstor/share/groupprojectWSXX/groupX/path/to/tool/bin
export PATH=${TOOL_HOME}:$PATH
```


## upload/download to/from server
You need to set up a ssh config on your own device. In your home directory,
there should be a directory ".ssh". In that directory if it's not already there,
you can create a file called "config" and write the following into it
(exchanging "your-username" with your username):

```
Host sshgw
Protocol 2
ControlPath ~/.ssh/socket-%r@%h:%p
hostname sshgw.cs.uni-tuebingen.de
user your-username

Host nog
Hostname nog.cs.uni-tuebingen.de
user your-username
ControlPath ~/.ssh/socket-%r@%h:%p
ProxyCommand ssh sshgw -W nog:22

Host collins
Hostname collins.cs.uni-tuebingen.de
user your-username
ControlPath ~/.ssh/socket-%r@%h:%p

ProxyCommand ssh sshgw -W collins:22
```

Then, you can simply enter the server using `ssh nog` / `ssh collins` (you still need to enter your password twice).
If you have this set up you can use scp like this:

```
scp path/to/file/you/want/to/transfer nog:/teachstor/../path/to/copy/to/
```

This also works the other way around if you want to copy something from the
server to your laptop. (Basically like "cp")



