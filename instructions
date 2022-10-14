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

Using git commit with the -a option will automatically run the git add command for files that are already being tracked. This allows the user to explicitly skip the git add step for files that are already tracked. It will still be necessary to use the git add command for files that are not yet being tracked.

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

