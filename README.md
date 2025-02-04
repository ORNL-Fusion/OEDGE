To clone the code:

You can access the OEDGE code from github using either the SSH or HTML interface. The SSH approach will require an appropriate ssh encryption key pair be set up and the public key uploaded to your account on github.

It is recommended to check out the repository into a directory called divimp since the scripts within the package expect the code to be in a directory called divimp. The default behaviour of git if "divimp" is not included on the command line is to place the repository into a directory called OEDGE. That will work fine but be aware that execution scripts within the repo will need to be modified for your local installation.

SSH:
git clone git@github.com:ORNL-Fusion/OEDGE.git divimp

OR

HTML:
git clone https://github.com/ORNL-Fusion/OEDGE.git divimp

This process will copy the repo and place you on the master branch. The remote OEDGE.git repository will be associated with the name "origin". To list all of the remote OEDGE branches use the following command:

git branch -r  or git branch --remotes

The typical output of this command will look like the following:
  origin/HEAD -> origin/master
  origin/aalto
  origin/asipp
  origin/asipp-zang
  origin/aug
  origin/cedric
  origin/diiid
  origin/diiid-greg
  origin/diiid-jerome

... and so on.


The command "git branch" --list or "git branch -l" by itself will list the branches available in your local repository.

If you already have an OEDGE branch which you are working on but not checked out locally you can establish a local branch to track the remote branch so that changes made locally can be pushed to the remote branch and so you can merge the changes of other collaborators to the remote branch into your local branch.

The following command creates a local remote tracking branch.

git checkout --track origin/<branch-name>

<branch-name> would be taken from the list of remote branches available in the repository. This will create a local branch which tracks the remote branch and will switch your current branch to the local copy of the remote branch. You can see which branch you are working on in your local repository using the git branch command. The branch with an asterix beside it is the current acive branch.

If collaborators wish to merge changes to the master branch then they will need to merge the master into their branch, fully test the branch for both their own workflows and a sample of workflows used by other users and then issue a pull request on github to merge these changes into the master branch. 


GIT reference:

A good reference for GIT is the Pro GIT book. This contains details of all the git commands. 

https://git-scm.com/book/en/v2

Some useful git cheat sheets can be found here:

Atlassian:
https://www.atlassian.com/git/tutorials/atlassian-git-cheatsheet
file:///home/david/Downloads/SWTM-2088_Atlassian-Git-Cheatsheet.pdf

Github:
https://education.github.com/git-cheat-sheet-education.pdf

Gitlab:
https://about.gitlab.com/images/press/git-cheat-sheet.pdf


A brief summary of useful git commands can be found below. 


Useful git commands:


1) git branch

The git branch command can also be used to create a new branch or track a remote branch. The following command creates a local branch with the name <branch name>. If <branch name> is on the remote repository then this will create a remote tracking branch from the repository (the remote repository will have the default alias of "origin").

git branch <branch name>

List local branches
git branch -l 

List remote branches
git branch -r

List all branches including local and remote
git branch -a

Create a new branch
git branch -b <branch name>

Delete a branch
git branch -d <branch name>


2) git switch /  git checkout <branch name> / git checkout -b <branch name>

You can change the branch you are working on locally using the git switch command. For example:

git switch <branch name>

The git checkout command is similar to git switch. It will change tracking to the specified branch if it exists. If you want to create the branch and change to it at the same time then the -b flag is needed.


3) git add <file name>

When changes are made to a file or a new file is added, the git add command is used to mark the file for updating. Using this command stages the files which will be committed on the next git commit command.

git add <file name>

This adds the file to the list of files to be updated in the next commit.

git add -f <filename> can be used to add files that are in the ignored list (.gitignore note at end) to be tracked by git.

4) git reset <file name>

This command will unstage a file that has been added with the git add command but it retains the changes made in the working directory copy.

5) git rm

This will remove a file from the local working directory and will stage the command to remove the file from repository tracking. The next time git commit is run the file will be removed from tracking in the repository. 

6) git diff

This command reports the differences between versions of the file.

"git diff" reports the difference of what is changed but not yet staged.
"git diff --staged" reports the differences of what is staged but not yet committed.
 
7) git commit -m "Description of commit"  / git commit -a

The git commit command will commit all the currently staged files. This effectively saves the changes to the local repository and logs the changes. The description is included in the update to summarize the changes. git commit without the -m option will open an editor to allow a note describing the commit to be created. 

Using git commit with the -a option will automatically run the git add command for files that are already being tracked and have been modified. This allows the user to explicitly skip the git add step for files that are already tracked. It will still be necessary to use the git add command for files that are not yet being tracked.

8) git status

The status of files in the local working copy can be found using this command. New files that are not yet tracked will be identified as untracked.  

9) git stash


Any files that have been staged but not committed changes will need to be stashed using the git stash command before switching branches or updating the branch with remote changes using the git pull or git fetch command. After the local branch has been updated, the changes that have been stashed can be re-applied using:

git stash pop

The following git stash commands can also be useful for managing these temporary changes.

List the stack order of stashed changes
git stash list

Write the saved changes from the stash to the local working copy branch
git stash pop

Discard the changes from the top of the stashed stack of changes.
git stash drop

10) git merge <branch name>

Merge the changes made in branch name into the current local branch. The merging process attempts to be automated but various change conflicts can occur between local and remote changes or changes between local branches that will require manual merge conflict resolution.

11) git pull / git fetch

The git fetch command will get changes made to the various branches on the remote repository and copy them to the local repository. However, these changes will not be applied to any of the local branches. The git pull command fetches changes from the matching remote branch AND merges those changes into the local branch. 

12) git push

This command pushes local changes to the branch to the remote tracking branch making them available to anyone else working on the same branch.

13) git checkout <-b branch> / git checkout <-- file>     (The <> surround optional arguments)

git checkout -- file

This command is used to revert changes in an uncommited tracked file in the local repo. It replaces the local file with an unmodified version from the repo. The changes that were in the local copy of the file are unrecoverable. This is useful for getting rid of debugging or other changes made locally that should not be kept. However, care should be taken when using this since the file is over written and all changes are lost. 

git checkout -b branch

This command changes the working directory to the specified branch. If the branch does not exist then git will create it. 

git checkout

This command actually does nothing :). From the git reference manual: 

"You could omit <branch>, in which case the command degenerates to "check out the current branch", which is a glorified no-op with rather expensive side-effects to show only the tracking information, if exists, for the current branch."


-----

Note: .gitignore

A .gitignore file has been added to the DIVIMP/OEDGE repo.

Most users utilize the DIVIMP source tree for code runs with input in the data directory and output in the results directory. The scripts in DIVIMP default to these locations though the scripts can be easily configured to have separate tree structures for the user files and the executable code.

However, since many users utilize the repo for code execution, a .gitignore file has been added so that git will ignore input and output files as well as the local libraries, executable files and various other buld files (see list below). Files in these locations that are already tracked will remain tracked. However, untracked files in these locations will not be shown in git status or accidentally added to the repo if a "git add ." or "git add -A" command is executed.

The current contents of the .gitignore file are:
*.o
*.mod
libsrc/**
local/
*~
results/**
data/**
cases/**
shots/**
div6/Makefile
out6/Makefile
eirene07/Makefile
eirene99/Makefile
rundiv
runout

Compiler files: .o and .mod
Library build files if built from source in libsrc: libsrc/**
Local library directories: local/
Emacs temporary files: *~
Content of the results, data, cases and shots directories: results/** data/** cases/** shots/**
Specific files for which local copies are created and configured: div6/Makefile out6/Makefile eirene07/Makefile eirene99/Makefile rundiv runout
Executable file names: div6/div6O* out6/out6O* eirene07/eirene eirene99/eirene triangle/showme triangle/triangle

In order to start tracking an ignored file, the -f or --force option must be used on the git add command:
git add -f <file name>
This should add the ignored file to the repo.


------

Notes: GIT idiosyncracies

1) GIT will not keep empty directories. It works only with files. If you want to create a directory tree within the repo to contain files generated locally (e.g. results) then it is suggested to create an empty file called '.gitkeep' in the otherwise empty directory and then add this file to the repo using "git add .gitkeep". GIT then tracks this file and will create the directory when cloned/checked out. 
