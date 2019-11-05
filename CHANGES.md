# Changes

## Next

* Changes required for use with GRCh38 compatible Barrenberg R
* **some reference files for for chromosome X will require renaming**
* `mv 1000_genomes_GC_corr_chr_23.txt.gz 1000_genomes_GC_corr_chr_X.txt.gz`
* `mv 1000genomesAlleles2012_chr23.txt 1000genomesAlleles2012_chrX.txt`
* `mv 1000genomesloci2012_chr23.txt 1000genomesloci2012_chrX.txt`

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
