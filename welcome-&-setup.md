# Welcome to the RNA-seq Data Analysis workshop #

Before you attend the workshop there are a couple of things we would like you to do to get setup so you are able to participate in all sections of the workshop.  

For those of you that indicated you did not have an account on the [Dartmouth Discovery cluster](https://rc.dartmouth.edu/index.php/discovery-overview/) you should have received an email explaining how to set that up, please make sure this is done and you are able to log into your account **BEFORE** the workshop begins. **YOU WILL NEED A DISCOVERY ACCOUNT!**


**Install VPN Client** 
You will also need to install Dartmouth College's VPN client to access Dartmouth compute resources. Once you have your netID you can navigate [here](https://services.dartmouth.edu/TDClient/1806/Portal/KB/ArticleDet?ID=72204) and follow the instructions to install the VPN client.

## How to do this workshop ##

We will be using a dataset downloaded from the [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra), a public repository of genomic data. This dataset comes from [this paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0099625), and was collected from human airway smooth muscle cells to test gene pathways effected by exposure to Glucocorticoid drugs, which have been historically used for their anti-inflammatory effects to treat asthma.

All the teaching materials are located within the [GitHub repository](https://github.com/Dartmouth-Data-Analytics-Core/RNA-seq_2023). We suggest you bookmark this page so you can easily get back to this repository each day.

Lessons are provided in Markdown format (files with extension  `(.md)`) and also contain *'code chunks'* that you will use to perform the analysis of this dataset. The majority of analysis will be performed using a terminal application or emulator, with an open `ssh` connection to the Discovery cluster. You may copy and paste the code chunks into your terminal window to perform analyses, or type them out by hand.  

If you wish to edit, modify or save the code in its own file (as you would in a real analysis), this can be achieved using a Text Editor application. Numerous free text editors exist, such as [Sublime Text](https://www.sublimetext.com/), [BBEdit](https://www.barebones.com/products/bbedit/), and [Atom](https://atom.io/). Experimenting with your code in a text editor is an excellent way to learn, as well as document your work on a project.

The terminal application you use will depend on the operating system you are using. If you do not already have a terminal emulator on your machine, please download one. Below are some recommendations for different operating systems.

Operating system| Terminal emulators
---|---
Mac| Terminal (comes pre-installed) <br> [iTerm2](https://iterm2.com/)
Windows| [MobaXterm](https://mobaxterm.mobatek.net/download.html) <br> [PuTTY](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html) <br> [iTerm2](https://iterm2.com/)
Linux| Konsole, Terminal, etc. (should be pre-installed but depends on the desktop environment you are running)


-------------

## The Dataset ##

For this workshop we will be using a published dataset as part of a study to determine the effects of Glucocorticoid drugs of human airway smooth muscle cells ([published in PLoS One](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0099625)). In the study, four cell lines were treated with either a control vehicle (untreated), dexamethasone (dex), albuterol (alb), or both dexamethasone and albuterol (co-treated) for 18 hours before transcriptomes were extracted.

The cartoon below provides an overview of the experimental design:

<p align="left">
<img src="figures/experiment.png" alt="overview"
	title="" width="60%" height="60%" />
</p>

The raw data was downloaded from the [Sequence Read Archive, SRA ](https://www.ncbi.nlm.nih.gov/sra?) a large NCBI database of high-throughput sequencing (HTS) data, and processed to generate a gene expression matrix of raw counts. These data are available under SRA accession [SRP033351](https://www.ncbi.nlm.nih.gov/sra?term=SRP033351).

Normalized data are also available from the [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/), another NCBI database, used to store processed HTS datasets, under accession [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778).


-------------

## Install the Integrative Genomics Viewer (IGV)

We will be using the [Integrative Genomics Viewer (IGV)](https://software.broadinstitute.org/software/igv/), a genome browser produced by researchers at the Broad Institute, to explore and visualize RNA-seq data.

<img src="figures/igv.png" height="100" width="100"/>

You will need to download and install the IGV desktop application for your operating system before the workshop begins. The latest versions of IGV can be found at their [downloads page](https://software.broadinstitute.org/software/igv/download). After installing IGV, try opening the application on your computer to confirm the installation was successful.

---
## Installing an SFTP client ##

This is optional, but for those of you that are new to the command line, an SFTP client might be an easier way to move files between the HPC environment and your local machine. SFTP stands for Secure File Transfer Protocol and will enable you to drag and drop files as you might in a finder window between your local machine and a remote location, rather than using the command line.

Several free SFTP clients exist, such as [FileZilla](https://filezilla-project.org/download.php?show_all=1) [WinSCP](https://winscp.net/eng/index.phpby), and [Cyberduck](https://cyberduck.io/), among others.

---


## Setting up a Conda Environment ##

Conda is an open source package and environment manager that runs on Windows, MacOS and Linux. Conda allows you to install and update software packages as well as organize them efficiently into environments that you can switch between to manage software collections and versions.

We will be using Conda to make sure everyone has the required software to perform the analyses included in the workshop. To start using Conda on Discovery, open your terminal application and start an ssh connection using your username & password:

```bash
# Establish the secure shell connection
#### REPLACE 'netID' WITH THE ID ASSOCIATED WITH YOUR DISCOVERY ACCOUNT
ssh netID@discovery7.dartmouth.edu

# Enter your password at the prompt (when you type no characters will show up to preserve privacy)
netID@discovery7.dartmouth.edu password:

# You're in!
(base) [netID@discovery7 ~]$
```

Then run the following command:

```bash
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
```

We recommend that you add the above line of code to your `.bashrc` file in your home directory, otherwise you will need to run this command each time you start a new session on discovery. To do this use the following code:

```bash
# navigate to your home directory
cd ~
# open the .bashrc file that is there
nano .bashrc
```
This will open the existing `.bashrc` file use the down arrow to move the cursor to the bottom of the file and paste `source /optnfs/common/miniconda3/etc/profile.d/conda.sh`. Then use the `ctrl + x` keys to exit the nano text editor, type `Y` to save the changes you made, and hit `return` to save the file to the same name (`.bashrc`).

Now run the following command to create a .conda/ directory in your home drive to store all of your personal conda environments. You only have to run this command once to make this directory, so it does not need to be added to your `.bashrc` file.

```bash
cd ~
mkdir -p .conda/pkgs/cache .conda/envs
```

When you are ready activate the conda environment for the workshop, use the following command:

```bash
conda activate /dartfs-hpc/scratch/rnaseq1/envs/rnaseq1
```

You will see that the activate command has worked when it reads `rnaseq1` rather than `base` to the left of the prompt.

When you are finished using a conda environment, it is good practice to deactivate your session with the following command.

```bash
conda deactivate
```

That's it! This conda environment contains all the software you will need during the workshop. If you run into issues with the setup, please reach out to us at *DataAnalyticsCore@groups.dartmouth.edu* and someone will be in touch to assist you.

**NOTE:** Dartmouth's Research Computing team also provides instructions for getting started with Conda on discovery, which you can find [here](https://services.dartmouth.edu/TDClient/1806/Portal/KB/ArticleDet?ID=72888).  


---
## Write a test file to the scratch space ##

We've created a folder on the scratch space for this workshop where everyone can write 'Hello' once they've completed these welcome and setup instructions.  To write your own file, use the following code, replacing 'xyz' and 'your_name' with your own.  The quotation marks and spaces are important!  This will create a record of how many of us have successfully logged in to Discovery and finished the welcome and setup tasks:  

```bash
echo "Hello from xyz" > /dartfs-hpc/scratch/rnaseq1/welcome/your_name.txt
```
So, for example, Tim would write:  
```bash
echo "Hello from Tim" > /dartfs-hpc/scratch/rnaseq1/welcome/tim_sullivan.txt
```

----

## Install R & RStudio ##

### R
For all the analysis in the workshop, we will be using R, a free open source programming language and statistical software environment used extensively in bioinformatics. R is also a powerful way to create high quality graphics.

Visit the *R-Project* homepage [here](https://www.r-project.org/) and download a R version (4.0.3 or greater) that is appropriate for your operating system.

<p align="left">
<img src="figures/r-logo.png" alt="overview"
	title="" width="25%" height="25%" />
</p>

### RStudio

To help use R efficiently, we will also be using [RStudio](https://rstudio.com/), an IDE (Integrated Development Environment) for R built to consolidate different aspects of writing, executing, and evaluating computer code. Without an IDE, these aspects of programming would need to be performed in different applications,  reducing productivity.  

<p align="left">
<img src="figures/r-studio-logo.png" alt="overview"
	title="" width="40%" height="40%" />
</p>

**Basic features of the RStudio IDE include:**  
- console for submitting code to
- syntax-highlighting editor used for writing R-scripts
- windows for environment management, data visualization, and debugging
- facilities for version control & project management

Navigate to the [RStudio website](https://rstudio.com/) and download the appropriate version for your operating system.

-------------

## Install required R-packages

Beyond the basic functionality included in R's standard distribution, an enormous number of packages designed to extend R's functionality for specific applications exists, representing one of R's core strengths.

Most R-packages are obtained from one of two package repositories:
- [CRAN (The Comprehensive R Network)](https://cran.r-project.org/)
- [Bioconductor](https://www.bioconductor.org/)

During the workshop we will be using a number of packages from both CRAN and Bioconductor. Once you have installed R and RStudio, open RStudio (or R) and copy & paste the following code chunk into the console. This will prompt R to download and install the specified packages. We have demonstrated this for you in [this video](https://youtu.be/UtZHS-q7buI).



```r
if (!any(rownames(installed.packages()) == "tximport")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("tximport")
}
library(tximport)

if (!any(rownames(installed.packages()) == "DESeq2")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("DESeq2")
}
library(DESeq2)

if (!any(rownames(installed.packages()) == "biomaRt")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("biomaRt")
}
library(biomaRt)

if (!any(rownames(installed.packages()) == "vsn")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("vsn")
}
library(vsn)

if (!any(rownames(installed.packages()) == "ComplexHeatmap")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("ComplexHeatmap")
}
library(ComplexHeatmap)

if (!any(rownames(installed.packages()) == "readr")){
      install.packages("readr")
}
library(readr)

if (!any(rownames(installed.packages()) == "ggrepel")){
      install.packages("ggrepel")
}
library(ggrepel)

if (!any(rownames(installed.packages()) == "rlang")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("rlang")
}
library(rlang)

if (!any(rownames(installed.packages()) == "EnhancedVolcano")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("EnhancedVolcano")
}
library(EnhancedVolcano)

if (!any(rownames(installed.packages()) == "apeglm")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("apeglm")
}
library(apeglm)

if (!any(rownames(installed.packages()) == "dplyr")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("dplyr")
}
library(dplyr)

if (!any(rownames(installed.packages()) == "ggplot2")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("ggplot2")
}
library(ggplot2)

if (!any(rownames(installed.packages()) == "pheatmap")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("pheatmap")
}
library(pheatmap)

if (!any(rownames(installed.packages()) == "gplots")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("gplots")
}
library(gplots)

if (!any(rownames(installed.packages()) == "RColorBrewer")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("RColorBrewer")
}
library(RColorBrewer)

if (!any(rownames(installed.packages()) == "circlize")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("circlize")
}
library(circlize)

sessionInfo()

```

----

#### Troubleshooting R package installation

If one of the packages fails to install correctly due to the absence of another package that is not installed by the above code chunk, simply install the missing package, then re-run the failed code from above.

- If the missing package is from the CRAN repository, use `install.packages(XXX)`.

- If the missing package is from Bioconductor, use `  BiocManager::install("XXX")`.

If you are unable to resolve the issue, please reach out to us for assistance **before** the workshop at DataAnalyticsCore@groups.dartmouth.edu.
