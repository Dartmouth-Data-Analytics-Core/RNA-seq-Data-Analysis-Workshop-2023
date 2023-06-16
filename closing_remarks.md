# Closing remarks

## Workshop goals:
- Understand the basic principles of a differential expression (DE) analysis using RNA-seq data
- Develop a working understanding of the fundamental statistics fro typical DE analysis methods
- Perform a DE analysis using R/Bioconductor packages
- Learn how to explore the results and make robust insights from your data

----------

## DE analysis overview
![](figures/day2_summary.png)

----------

## Final take-aways from the workshop
- Spend the time to plan, consult, practice, (and money) to generate high quality data sets

- If you are going to do a lot of Bioinformatics, you should get really good at the command-line (Bash), otherwise, pre-processing will be slow & painful

- Identify, understand, and check key QC metrics to ensure the quality of your results

- Spend the time to plan your experiment well (enough replicates, appropriate sequencing depth, etc.). No analysis can rescue a bas dataset.

- If you will perform differential expression analysis regularly, you should build your experience in R, as well as you knowledge of fundamental statistical concepts

- Make sure you understand what the core functions (e.g. `DESeq2`) are doing by consulting the original manuscripts documentation. Failure to understand what these functions do may lead you to perform your analysis incorrectly.

- Correct for multiple testing!

----------

## How to consolidate your learning:
- Revisit the workshop materials a week or two after the workshop, and re-run the analysis code from scratch

- Edit/annotate the code, run sub-sections, and read the `help` pages for important functions

- Practice with the complete dataset (all chromosomes), that is available to you for approx. 1 month on discovery in /scratch/rnaseq1/data/. This will give you experience running data from the entire genome, and an appreciation for the computational resources and required time to complete these tasks.

- Read the methods sections of published papers that perform differential expression analyses. This will help you to understand how many of the concepts we have discussed are applied and reported in practice

- Read reviews like [this one](https://pubmed.ncbi.nlm.nih.gov/31341269/) from Stark *et al*, 2019, *Nat. Rev. Genetics*, `RNA Sequencing: The Teenage Years`.

- Ask us questions! (Bioinformatics office hours: https://dartmouth.zoom.us/s/96998379866, Friday's at 1-2 pm, password: *bioinfo*)

----------

## High performance computing (HPC)


During the workshop, we all logged onto Discovery and ran our workflow interactively. Generally, this type of 'interactive' work is not encouraged on Discovery, and is better performed using other servers such as Andes & polaris (see article [here](https://rc.dartmouth.edu/index.php/hrf_faq/how-do-i-run-interactive-jobs/) on this topic from research computing).

However, working interactively with genomics data can be quite slow since many operations will need to be run in sequence across multiple samples. As an alternative, you can use the scheduler system on discovery, that controls how jobs are submitted and managed on the Discovery cluster.

Using the scheduler will allow you to make use of specific resources (such as high memory nodes etc.) and streamline your workflows more efficiently. Dartmouth just transitioned to using [Slurm](https://services.dartmouth.edu/TDClient/1806/Portal/KB/ArticleDet?ID=132625).

We encourage you to get comfortable using the scheduler and submitting jobs using an HPC system. [Research Computing](https://rc.dartmouth.edu/) has a lot of great material regarding using Discovery on their website.

----------

## Suggested reading:

Reading manuscripts that use RNA-seq, or reviews specifically focused on RNA-seq are excellent ways to further consolidate your learning.

In addition, reading the original manuscripts behind common tools will improve your understanding of how that tool works, and allow you to leverage more complicated options and implementations of that tool when required.

Below we provide some suggested reading to help get you on your way:

#### Review articles
- [Stark *et al*, 2019, *Nat. Rev. Genetics*.](https://pubmed.ncbi.nlm.nih.gov/31341269/) `RNA Sequencing: The Teenage Years`
- [Conesa *et al*, 2016, *Genome Biology*.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8) `A survey of best practices for RNA-seq data analysis`
- [Wang, *et al*, 2009, *Nat. Rev. Genetics*.](https://www.nature.com/articles/nrg2484) `RNA-Seq: a revolutionary tool for transcriptomics`
- [Cresko Lab, Univeristy of Oregon.](https://rnaseq.uoregon.edu/) `RNA-seqlopedia
: provides an overview of RNA-seq and of the choices necessary to carry out a successful RNA-seq experiment.`

#### Original manuscripts: Popular RNA-seq tools
- [Cutadapt:](http://journal.embnet.org/index.php/embnetjournal/article/view/200) `Cutadapt Removes Adapter Sequences From High-Throughput Sequencing Reads`.
- [STAR:](https://academic.oup.com/bioinformatics/article/29/1/15/272537) `STAR: ultrafast universal RNA-seq aligner`
- [HISAT2:](https://www.nature.com/articles/s41587-019-0201-4) `Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype`
- [Bowtie2:](https://www.nature.com/articles/nmeth.1923)`Fast gapped-read alignment with Bowtie 2`
- [HTSeq-count:](https://academic.oup.com/bioinformatics/article/31/2/166/2366196)`HTSeq—a Python framework to work with high-throughput sequencing data `
- [DESeq2:](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)`Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2`

----------

## Post DE analysis

After completing a DE analysis, we are usually left with a handful of genes that we wish to extract further meaning from. Depending on our hypothesis and what data is available to us, there are several ways to do this.

### 1. Integrative genomics

If you also collected other genomic data of a different modality ion the same samples (e.g. WGS-/WES, ChIP-seq, ATAC-seq), or an appropriate public dataset exists, you may choose to integrate your DE results with these data. This approach is referred to as *integrative genomics* and allows you to make insights that you would be unable to with a single data type.

For example, if you collected ChIP-seq for a specific transcription factor (TF) with paired RNA-seq data, you may wish to use your significant DEGs to identify genes whose expression is turned on/off by this TF under a treatment condition.

<p align="center">
  <img src="figures/int-gen.png" height="80%" width="80%"/>
</p>

### 2. Gene ontology (GO) & pathway analyses
Unless very few significant DEGs were identified in your analysis, it can be difficult to extract biological insights from a long list of genes. GO and pathway analysis methods represent a diverse collection of approaches that seek to prioritize sets of functionally related genes that are enriched in a list of DEGs.

<p align="center">
  <img src="figures/go-methods.png" height="80%" width="80%"/>
</p>

Many tools and methodologies exist for performing GO & pathway enrichment analysis. These tools make use of varied statistical approaches, some of which were designed for specific applications (such as for microarray data, not RNA-seq, e.g. GSEA), therefore selecting the most appropriate tool for your analysis can be non-trivial. We encourage you to read more about these methods and review your analysis plan with an expert if you plan to use these in you research.

Some suggested reading regarding gene ontology and pathway analysis approaches:  
- [Gene set analysis approaches for RNA-seq data: performance evaluation and application guideline. *Briefings in Bioinformatics.* 2016.](https://doi.org/10.1093/bib/bbv069)

- [Ten Years of Pathway Analysis: Current Approaches and Outstanding Challenges. *PLoS Computational Biology.* 2012.](https://doi.org/10.1371/journal.pcbi.1002375)

- [Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS* 2005.](https://doi.org/10.1073/pnas.0506580102) (the original GSEA paper)

- [Gene Set Enrichment Analysis Made Simple. *Stat Methods Med Red.* 2009.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3134237/)

-----------

## Feedback

We ask that you complete the survey that will be sent out over email so that we can gauge what worked well, and what we need to improve for the future. If you have additional thoughts that were not addressed in the survey, please feel free to contact any one of us, or reach out to the DAC email directly (*DataAnalyticsCore@groups.dartmouth.edu*).

<img src="figures/logo.jpg" width="250" height="140" >

This workshop will be offered again, in addition to our other bioinformatics workshop offerings (details will be posted on [CQB website](https://sites.dartmouth.edu/cqb/)). If you have suggestions for workshops you would like to see, please let us know!

---------

## Bioinformatics consultations

Please reach out to us with questions related to content from this workshop, or for analysis consultations. We also host office hours and consultations by request.

