# To-do
 - Acquire data
 - 

# Code
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