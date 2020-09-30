# Changes

## 3.7.0
 * Updated R Battenberg code to v2.2.9

## 3.6.0
 * Added option the reference download script to use a pre-downloaded archive
 * Updated Dockerfile base image
 * Removed redundant logic from setup script

## 3.5.3

* Fix parsing of subclones file when generating final outputs
* This release points to cancerit fork of core battenberg algorithm [here](https://github.com/cancerit/battenberg/feature/grch38)
  * v2.3.1 - faster PCF function

## 3.5.2

* Fix parsing of subclones file during refit

## 3.5.1

* Fix #107, again... solves the problem of Rprofile not being found when container invoked by tools that mess with $HOME

## 3.5.0

* Changes required for use with GRCh38 compatible Battenberg R
* **some reference files for for chromosome X will require renaming in GRCh37**
* `mv 1000_genomes_GC_corr_chr_23.txt.gz 1000_genomes_GC_corr_chr_X.txt.gz`
* `mv 1000genomesAlleles2012_chr23.txt 1000genomesAlleles2012_chrX.txt`
* `mv 1000genomesloci2012_chr23.txt 1000genomesloci2012_chrX.txt`
* This release points to cancerit fork of core battenberg algorithm [here](https://github.com/cancerit/battenberg/feature/grch38)

## 3.4.0

* Placed into a container.

## 3.3.1

* Fixed bug when creating split loci files which would occasionally count wrongly due to rounding errors. Fixed the cause of the rounding error and added a sort to the hash.
* Added check that the number of split loci files is at least the number of required contigs

## 3.3.0

* Upgrade to Battenberg [v2.2.8](https://github.com/Wedge-Oxford/battenberg/releases/tag/v2.2.8).
* Upgrade to ASCAT [2.5.1](https://github.com/Crick-CancerGenomics/ascat/releases/tag/v2.5.1).
* Docs now indicate you should run the R install scripts.
* Email and licence updates

## 3.2.2

* Fix to README rendering of link

## 3.2.1

* Updates to documentation

## 3.2.0

* Install Battenberg R code from GitHub

## 3.1.0

* Use newer alleleCount that includes faster 'dense snps' option
* Minor modification to README.md layout

## 3.0.2

* Moved call to create 1000genome loci filename hash to before use of threads to avoid glob bug in perl 5.16.3

## 3.0.1

* Bug fix when there are regions which have no copynumber values set which caused bb_vcf_to_ascat_cn.pl to crash

## 3.0.0

* Use verion 2.2.5 of the battenberg R package
* Removed R files from repository
