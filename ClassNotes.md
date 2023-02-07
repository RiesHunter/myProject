# Class Notes
## 230126
### Created [GitHub repository](https://github.com/RiesHunter/myProject). 
Updated git, configured git config, and cloned the class and personal repos:
```shell
brew upgrade git
git config --global user.name 'RiesHunter'
git config --global user.email 'Hries@wisc.edu'
git clone https://github.com/crsl4/phylogenetics-class
git clone https://github.com/RiesHunter/myProject
```

### Updated README.md to test version control
added: "Test update"
```shell
git add .
git commit -m "updated readme"
git push
```

### Created test directories and placed foo.txt and bar.txt in each directory with:
```shell
mkdir scripts data results figures manuscript
find . -type d -exec touch {}/foo.txt \;
find . -type d -exec touch {}/bar.txt \;
rm .git/foo.txt .git/bar.txt
git add .
git commit -m "test directories"
git push
```
### Created Repro.md file and pushed
Installed Copy Markdown as HTML:Name: 
    Copy Markdown as HTML
    Id: jerriepelser.copy-markdown-as-html
    Description: Copies the selected text from a markdown document to the clipboard as HTML.
    Version: 1.1.0
    Publisher: Jerrie Pelser
    VS Marketplace Link: https://marketplace.visualstudio.com/items?itemName=jerriepelser.copy-markdown-as-html)

Used [this website](https://www.freecodecamp.org/news/how-to-open-visual-studio-code-from-your-terminal/) to enable `code` command for Terminal. I also changed the default application for README.md and Repro.md to `VS Code`, so I can use `open` for these files specifically.

```shell
touch Repro.md
code Repro.md # added lots here
git add .; git commit -m "Added Repro.md file"; git push
```
### More reformatting and a bit of cleaning
```shell
#two pushes
git add .; git commit -m "Reformatted Repro.md"; git push
rm foo.txt bar.txt
git add .; git commit -m "More reformatting"; git push
```

## 230131
### Cloning course repo
```shell
#clicked "fork" [here](https://github.com/crsl4/phylogenetics-class), creating my own fork of the repo [here](https://github.com/RiesHunter/phylogenetics-class). This was nearly instantaneously completed
#this allows me to pull her new commits to my fork
#I had to remove my old copy of /phylogenetics-class that I had cloned from her repo
cd /Users/rieshunter/Documents/bioinformatics
git clone https://github.com/RiesHunter/phylogenetics-class #This will take a long while to clone a local version
cd phylogenetics-class
git remote -v # will tell us where we are in git
git remote add upstream https://github.com/crsl4/phylogenetics-class.git #connects my local clone and her origin
#output:
    #origin	https://github.com/RiesHunter/phylogenetics-class (fetch)
    #origin	https://github.com/RiesHunter/phylogenetics-class (push)
    #upstream	https://github.com/crsl4/phylogenetics-class.git (fetch)
    #upstream	https://github.com/crsl4/phylogenetics-class.git (push)
```
### Course repo is massive, so we're redo-ing it with a smaller repo
```shell
# fork [crsl4/phylo-class-social](https://github.com/crsl4/phylo-class-social)

# clone my fork to local
cd /Users/rieshunter/Documents/bioinformatics
git clone https://github.com/RiesHunter/phylo-class-social.git
cd phylo-class-social
git remote -v #clone not linked to origin

# link fork to origin
git remote add upstream https://github.com/crsl4/phylo-class-social.git #"upstream" is typical, but not a required name—you can name it whatever you'd like!
git remote -v #clone now linked to origin

# pulling upstream changes
# exploring potential errors
## make sure to pull from repo before pushing
```

### In-class exercise
```shell
# change best-books.md file
## - Hunter: Invisible Man by Ralph Ellison
git config pull.rebase false

git add .; git commit -mv "Hunter removed an empty directory"
git pull upstream master -v
git push -v
# have to manually do a pull request in GitHub
```

### Commit Repro.md for 230131
```shell
git add . -v; git commit -m "Formatting update: Updated repro from fork day"; git push -v
```

## 230202
### To-do
- [Learn@Home: Why learn phylogenomics?](https://github.com/crsl4/phylogenetics-class/tree/master/lecture-notes/lecture2-learn-home.md)
- [Learn@Home: Sequencing](https://github.com/crsl4/phylogenetics-class/blob/master/lecture-notes/lecture4-learn-home.md)
- HW 2.1 (reading)
- Git [here](https://github.com/crsl4/phylogenetics-class/blob/master/exercises/hw-git.md)

### Reproducibility pop quiz
```shell
# added file Hunter.md with some text
git add .; git commit -m "Hunter added a file."; git push

## looks like you can reset your master branch to the upstream with the following command
git reset --hard upstream/master
# I ran this and then did the following code:
nano Hunter.md #put some text in there
git add .; git commit -m "Should be mergeable"; git push
# looks like pull request updated with most recent commit after this!
```

### Sequencing
#### MSA
MSA helps resolve the historical substitution, insertion, and deletion evolutionary events
- Very computationaly intensive

Exercise 1:
C-->G sub
    ACAT
    AGAT
C deletion; G insertion
    AC-AT
    A-GAT

Exercise 2:
What would be the alignment of the sequence
    ACATTA
if it evolves into
    TACA
and we know the following happens:
- deletion of the first two nucleotides AC
- Deletion of the second T
- substitution of T into C
- inversion of T at the front
Workspace:
    -ACATTA
    T--A-CA
without knowing the true evolutionary events, we would have created the alignment:
    -ACATTA
    TACA---

Defining the cost of each event: deletion, insertion, substitution
- Important for weighting what you value
- E.g., an indel would likely be really bad in an intron

Exercise 3:
How would you align `AACT` and `CTGG`?
Cost of gap = 1
Cost of sub = 3
Workspace: Cost = 1(4) + 3(0) = 4
    AACT--
    --CTGG
I chose this because subs are costly (3) and indels are cheap (1)

Exercise 4:
How would you align `AACT` and `CTGG`?
Cost of gap = 4
Cost of sub = 1
Workspace: Cost = 4(0) + 1(0) = 0
    AACT
    CTGG
I chose this because any single indel would cost 4; whereas, all sub = 4
Since a single indel can't make the number of subs 0, we won't do an indel and instead will just do all subs

Needleman-Wunsch algorithm
- Pairwise sequence alignment
- Smart algorithm that works recursively in smaller chunks
    - Basically pairwise recursive comparison scores

## 230207
### To-do
- [Needleman-Wunsch HW](https://github.com/crsl4/phylogenetics-class/blob/master/exercises/hw-needleman.md)

Changed "Repro.md" to "notebook-log.md" per class instructions
Updated .gitignore to ignore "ClassNotes.md"
Quick push before class: "Notebook name change and placed ClassNotes.md in .gitignore"
Second push before class: "Removed ClassNotes.md from .gitignore"