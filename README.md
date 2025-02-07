
![DNA_CoreLab_logo](https://github.com/user-attachments/assets/dc990a43-7b10-466b-bd4b-040308990ba7)

The primary purpose of this app is to create a user friendly way to take advantage of R packages for parentage assignment and population evaluation using genomic data.

## Installation:
To use this locally, install and load the following R packages.
**R Packages:**
```R
install.packages(c("shiny", "sequoia", "Rcpp", "kinship2",
                   "ggplot2", "vcfR", "poppr", "ape", "RColorBrewer",
                   "reshape2", "adegenet", "cowplot", "Cairo",
                   "shinyWidgets", "grDevices", "shinyjs", "shinythemes",
                   "markdown", "DT"))

#define vector of packages to load
some_packages <- c("shiny", "sequoia", "Rcpp", "kinship2", 
                   "ggplot2", "vcfR", "poppr", "ape", "RColorBrewer", 
                   "reshape2", "adegenet", "cowplot", "Cairo", 
                   "shinyWidgets", "grDevices", "shinyjs", "shinythemes",
                   "markdown", "DT")

#load all packages at once
lapply(some_packages, library, character.only=TRUE)
```

Then you can run the app using:
```R
runGitHub("ParentageApp", "samstroupe")
```

## Parentage Assignment:
The parentage assignment function of this app uses the R package [Sequoia](https://jiscah.github.io/index.html), also available on [Github](https://github.com/JiscaH/sequoia), to determine the best sire and dam for each offspring [(Huisman 2017)(https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12665). The data input is a genotype file in `plink.raw` format and a csv file with at least ID, Sex, and Birth Year of each sample, see below. [VCFtools](https://vcftools.github.io/man_latest.html) can be used to filter the dataset and convert to PLINK PED format, then [Plink](https://www.cog-genomics.org/plink2/) can be used to convert to the PLINK.raw format for this app. Here is an example of how to convert data on the command line:
```bash
# Convert VCF to PLINK PED
vcftools –vcf <input.vcf> –plink –out <plink_out>

# Convert PLINK PED to PLINK RAW
plink –file <plink_out> –recode A –out <plink_raw_out>
```

The Life History file is a csv file with an ID, Sex, BY (Birth Year) Column. Optionally, you can include three additional columns BY.min, BY.max, and Year.last.

- ID: Sample identifier.
- Sex: 1 = female, 2 = male, 3=unknown, 4=hermaphrodites. All other numbers, letters, or NA = unknown
- BirthYear: Year of birth/hatching/germination. NA if birth year is unknown.
- BY.min: minimum birth year, only used if BirthYear is missing.
- BY.max: maximum birth year, only used if BirthYear is missing.
- Year.last: Last year in which individual could have had offspring. E.g. the year before death for females, the year after death for males, or year an animal was sold/transferred.

| ID | Sex | BY | BY.min | BY.max | Year.last |
| ----------- | ----------- | ----------- | ----------- | ----------- | ----------- |
| Sample_01 | 2 | 2005 | |||
| Sample_02 | 1 | 2008 | || 2022|
| Sample_03 | 1 | 2017 | |||
| Sample_04 | NA | 2020 | |||
| Sample_05 | 2 | NA |2014|2016||
| Sample_06 | 1 | 2020 ||||

## Population Evaluation:
Various packages installed above were used to evaluate population genetics data with this app. The purpose of the population evaluation functions are to quickly and easily identify variation and population substructure among samples. Many of the packages used were developed by the Grünwald lab. They have a helpful tutorial [here](https://grunwaldlab.github.io/Population_Genetics_in_R/index.html) for explanation on the methods used in this app and more. The input for these packages is a filtered vcf with only biallelic SNPs. [VCFtools](https://vcftools.github.io/man_latest.html), [Plink](https://www.cog-genomics.org/plink2/), [BCFtools](https://samtools.github.io/bcftools/bcftools.html), and [GATK](https://gatk.broadinstitute.org/hc/en-us/categories/360002369672-Tool-Index) are all good sources on how to filter and manipulate genetic/genomic data. 

The population input file is a simple csv file with a column for ID and Population. This App only takes the two needed columns, but similar to the Life History file, additional columns may be present in this input as well. It is important to note the ‘ID’ and ‘Population’ column headers must be identical to the example below.

| ID | Population | Additional |
| ----------- | ----------- | ----------- |
| Sample_01 | Pop_01 | columns |
| Sample_02 | Pop_01 | may |
| Sample_03 | Pop_02 | be |
| Sample_04 | Pop_02 | present |
| Sample_05 | Pop_02 | in |
| Sample_06 | Pop_01 | file |

## References:
- Huisman, J. (2017). Pedigree reconstruction from SNP data: parentage assignment, sibship clustering and beyond. Molecular ecology resources, 17(5), 1009-1024.
- Niklaus J. Grünwald, Zhian N. Kamvar, Sydney E. Everhart, Javier F. Tabima, and Brian J. Knaus © 2017, Corvallis, Oregon, USA
