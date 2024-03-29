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

Used [this website](https://www.freecodecamp.org/news/how-to-open-visual-studio-code-from-your-terminal/) to enable `code` 
command for Terminal. I also changed the default application for README.md and Repro.md to `VS Code`, so I can use `open` 
for these files specifically.

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
#clicked "fork" [here](https://github.com/crsl4/phylogenetics-class), creating my own fork of the repo 
[here](https://github.com/RiesHunter/phylogenetics-class). This was nearly instantaneously completed
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
git remote add upstream https://github.com/crsl4/phylo-class-social.git #"upstream" is typical, but not a required name—you 
can name it whatever you'd like!
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
- [Learn@Home: Why learn 
phylogenomics?](https://github.com/crsl4/phylogenetics-class/tree/master/lecture-notes/lecture2-learn-home.md)
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

### Alignments continued
 - Costs of alignment stack as you work across the matrix
`F(i,j) = min{F()+cost(),F()+1,F()+1}`
  - F(i,j), where F_i is row i and F_j is column j

```shell
# Example: 
#  - Cost_sub = 1
#  - Cost_gap = 2
#           f0 f1 f2 f3 f4
#              b1 b2 b3 b4
#           _  A  T  C  G
#  f0     _ 0  2  4  6  8
#  f1 a1  T 2  1  2  4  
#  f2 a2  C 4            
#  f3 a3  A 6          
 
# Choose the lowest cost for a1b1:
# cost(a1,b1) +  F(0,0)           = 1 + 0 [*]
# cost(a1,b1) + (F(0,1) + F(0,0)) = 1 + 2 + 2
# cost(a1,b1) + (F(1,0) + F(0,0)) = 1 + 2 + 2

# Choose the lowest cost for a1b2:
# cost(a1,b2) +  F(0,1)           = 0 + 2 [*]
# cost(a1,__) + (F(0,2) + F(0,0)) = 2 + 2 + 2
# cost(__,b2) +  F(1,1)           = 2 + 1

# Choose the lowest cost for a1b3:
# cost(a1,b3) +  F(0,2)           = 1 + 4
# cost(a1,__) +  F(0,3)           = 2 + 6 + 2
# cost(__,b3) +  F(1,2)           = 2 + 2 [*]
```
Trace back arrows to construct the alignment: 
 - diag arrow = match a and b
 - Right arrow = b matches gap
 - Down arrow = a matches gap
  Solution: 
   `ATCG`
   `-TCA`

Other methods:
 - Progressive (ClustalW)
    - Most widely used
 - Consistency-based scoring (T-coffee)
    - Improvement over progressive by using a more strict score function
 - Iterative refinement algorithm (muscle)
    - Improvement over progressive by doing sequential alignments until convergence of score

Progressive alignment:
 - Compute rooted binary tree (guide tree) from pairwise distances
    - Can also be root of problems with bad guide trees (relies VERY heavily on guide tree)
    - Run Needleman-Wunsch for all sequences into a pariwise distance matrix
        - Very time consuming
    - Distance-to-tree function to create guide tree
 - Build MSA from the bottom (leaves) up (root)
        - Mis-alignments in leaves propagate to root
How do you align alignments?
 - "Profile" = frequency of each nucleotide at each site
    - Ignore gaps in frequency calculations
    - (a1,A) = frequency of A in first position for alignment a
    - (b1,A) = frequency of A in first position for alignment b
    - Use Needleman-Wunsch to align alignment matrix a to alignment matrix b
        `cost(a_i,b_j) = sum P(x of a_i) P(y of b_j)`
        - x will sequentially equal A, C, T, or G (same with y), but x≠y
        - cost will be from matrix value, not nucleotide (because we don't have nucs)

## 230214
Going to use JT McCrone's "Stochastic" dataset:
NCBI BioProject (accession no: PRJNA412631)

Installed sra-toolkit, clustalw, t-coffee
```shell
sudo apt install sra-toolkit
sudo apt install clustalw; clustalw # exit with "X"
sudo apt install t-coffee -y; t_coffee
sudo apt install muscle

# I created test.sh to test these packages
# ran the following for script permissions:
cd /home/rieshunter/GitHub/myProject/data/tool_test
chmod 755 test.sh
## added directory to PATH in .bashrc from $HOME
PATH="$PATH:/home/rieshunter/GitHub/myProject/data/tool_test"
```

full test.sh script:
```shell
#!/bin/bash

cwd=$(pwd)
data_dir="/home/rieshunter/GitHub/phylogenetics-class/data"
echo $data_dir
ls -lh $data_dir

## clustalw
clustalw -ALIGN -INFILE=${data_dir}/primatesAA.fasta -OUTFILE=${cwd}/primatesAA-aligned-clustalw.fasta -OUTPUT=FASTA

## t-coffee
t_coffee ${data_dir}/primatesAA.fasta
mv ${cwd}/primatesAA.aln ${cwd}/primatesAA-aligned-tcoffee_MSA.aln
mv ${cwd}/primatesAA.dnd ${cwd}/primatesAA-aligned-tcoffee_GUIDE_TREE.dnd
mv ${cwd}/primatesAA.html ${cwd}/primatesAA-aligned-tcoffee_MSA.html

## muscle
muscle -in ${data_dir}/primatesAA.fasta -out primatesAA-aligned-muscle.fasta
```


## 230221
Distance-based methods
 - Not character-based
 - Based on an explicit model of evolution
 - Create dissimilarity matrix between sequences ("p-distances")
 - Covert to evolutionary distance by correcting for multiple events per site
    - This is where the model of substitution comes in
    - Jukes & Cantor (1969): d_AB = -3/4 ln(1 - 4/3 * p_distance)
 - Create a distance tree
    - Don't need to search space of trees
        - Big advantage compared to parsimony or likelihood methods
        - Faster
    - Disadvantages
        - Approximate tree (may not be optimum; disadvantage of all inference methods)
    - Cluster analysis
        - "ultrametric trees" where end nodes are equidistant from root
        - Assumes a molecular clock
        - Uses WPGMA (simple average from ancestor) and UPGMA (averaged by n taxa in cluster)
        - Just a tree, no confidence interval
        - "Evolutionary distances" unit = substitutions
        - Might not take into account indels
    - Minimum evolution
        - Additive trees
            - Better for non-clock-like behavior
        - Algorithm = neighbor-joining
            - Input = matrix of pairwise distances
            - Sum of distances for every end node
            - Rate-correct distance matrix

## 230223
Parsimony-based methods
 - Does not rely on models of evolution
 - The tree that minimizes the amount of evolutionary change required to explain the data
 - Assumptions
    - Most effective when rate of evolution is slow (not always the assumption)
    - Can perform well under high rates of evolution as long as there are no pathological inequalities
    - Main assumption: independence among characters
 - Methodology
    - Determine the amount of character change required to explain the data
    - Search over all possible ultrametric tree topologies
        - Newick format (AKA parenthetical format)
        - Compute lowest parsimony score ancestral tree per site
        - Uses dynamic programming (breaking up complex task into manner smaller tasks)
    - Fitch algorithm
        - Faster solution
        - Random root
        - Assigns costs to mismatching nodes
            - Match = R + L
            - No match = R + L + 1
            - Tips of leaves equal zero initially
            (((C,C),T),G)
            (C,C)           = {C}       = 0 + 0 = 0
            ((C,C),T)       = {C,T}     = 0 + 0 + 1 = 1
            (((C,C),T),G)   = {C,T,G} = 1 + 0 + 1 = 2
 - "Long branch attraction" (Felsenstein zone) is the kiss of death for this method
    - Avoid
    - Long branches in additive trees will be grouped together in ultrametric trees

## 230302 
Models of evolution
 - Probability
    - Count number of observations for length of time
        - Have to account for variation!
        - Ballpark estimate
 - Different models have different assumptions:
    - Symmetric around the mean?
    - Variation around the mean?
 - Poisson model
    - P(X=k) = (lambda^k * e^-lambda) / k!
    - k = number of observations
        - "What is the probability of seeing k number of observations?"
    - Lambda = average
    - E.g.,
        - Average number of observations = 8.4
        - What is the porbability of seeing 5?
        - (8.4^5 * e^-8.4) / 5!
        - (41821.2 * 0.0002) / 120 
        - 8.4 / 120 = .078
    - Imperfect model
        - Single parameter (mean)
        - Assumes same model regardless of events
 - Assume:
    - The mutation process is the same at every branch of the tree
    - We assume sites evolve independently
        - We can focus on the mutation process between two sites
    - All sites evolve the same
        - We can choose any site that fits the model
 - Probability model for the evolution of an ancestral site into a descendant site
    - Number of mutations on time t is assumed to follow a Poisson distribution
    - P(X = k) = ( (µt)^k * e^(-µt) ) / k!
    - µ = rate of mutation events
 - Continuos-time Markov chain
    - Probabilities for the next event do not depend on the time point where the chain is
    - Continuous numerical time, not integer
    - Non-discrete
    - Independent of place in previous step
        - Forgets about the past
    - Probabilities of the next state only depend on the current state
    - Must choose substitution model Q
    - Math:
        - Ergodic
            - ?
        - Time-reversible
            - Same probability there and back
Models of substituion
 - Jukes-Cantor Model (JC69)
    - Setting all site-frequencies to be 1/4
    - Quick but likely less accurate
 - Felsenstein model (F81)
    - Does not assume observed frequencies of substitutinos
    - Empirically determine the four site-frequency paramters
    - Assumes rate that a nucleotide is observed is the rate of mutation to that nucleotide
 - General Time Reversible (GTR)
    - Very flexible and general, but complicated
    - Many parameters
    - Empirically determine the ten site-frequency paramters
        - 6 constants, 4 πs
    - Also uses constants for each site
