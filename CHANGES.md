### 3.0.2
* Moved call to create 1000genome loci filename hash to before use of threads to avoid glob bug in perl 5.16.3

### 3.0.1
* Bug fix when there are regions which have no copynumber values set which caused bb_vcf_to_ascat_cn.pl to crash

### 3.0.0
* Use verion 2.2.5 of the battenberg R package
* Removed R files from repository
