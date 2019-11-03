# cgpBattenberg

An installation helper, perl wrapper and the R program Battenberg which detects subclonality and
copy number in matched NGS data.

[![Quay Badge][quay-status]][quay-repo]

| Master                                        | Develop                                         |
| --------------------------------------------- | ----------------------------------------------- |
| [![Master Badge][travis-master]][travis-base] | [![Develop Badge][travis-develop]][travis-base] |

**_This is only suitable for WGS analysis_**.

<!-- TOC depthFrom:2 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Battenberg R code](#battenberg-r-code)
- [Docker, Singularity and Dockstore](#docker-singularity-and-dockstore)
- [Installation](#installation)
	- [Prerequisites](#prerequisites)
- [Program Run Instructions](#program-run-instructions)

<!-- /TOC -->

## Battenberg R code

The Battenberg R code is maintained in a separate repository [Wedge-Oxford/battenberg][bb-repo]
and this is where any questions or issues specific to the R code should be directed.

## Docker, Singularity and Dockstore

There is a pre-built image containing this codebase on quay.io.

* [dockstore-cgpwgs][ds-cgpwgs-git]: Contains additional tools for WGS analysis.

This was primarily designed for use with dockstore.org but can be used as normal containers.

The docker images are know to work correctly after import into a singularity image.

## Installation

The battenberg R files are installed automatically from the Battenberg GitHub repository found
[here][bb-repo]. The linked version is currently [`v2.2.5`][bb-ver-link].

Please install the following first:

1. [PCAP-core v2.1.3+][pcap-core-rel]
1. [alleleCount v3.3.1+][allele-count-rel]
1. [cgpVcf v2.0.1+][cgpvcf-rel]

Then execute:

```
setup.sh <install_to_folder> [X/lib/perl:Y/lib/perl]
cd Rsupport
./setupR.sh <install_to_folder>/R-libs
```

All of the items listed here use the same install method.

### Prerequisites

* Impute2 executables can be found [here][impute-exe]
  * Any impute related data for download
* BWA Mapped, indexed, duplicate marked/removed bam files, for both a matched normal and tumour sample
* Reference.fasta and index
* A file containing a list of contigs in the reference .fai to ignore

Some required data files are not included in the distribution but a script is included to generate these for you:

* Directory containing the 1000 genomes allele and loci data:
  * Generated using the included script ``download_generate_bberg_ref_files.pl``
* Impute info file ``impute_info.txt``
  * Generated using the included script ``download_generate_bberg_ref_files.pl``
* Prob loci file probloci.txt
  * Included: ``files/probloci.txt.gz``

Additionally, the wgs_gc_correction_1000g files need to be downloaded. These can be obtained from the Battenberg R code site [here][bb-ref].

## Program Run Instructions

For the most up to date usage instructions for the wrapper code please see the command line help:

    battenberg.pl -h

Please check the [wiki][cgpbb-wiki] for common problems before raising any issues.

# LICENCE

```
Copyright (c) 2014-2018 Genome Research Ltd.

Author: Cancer Genome Project <cgphelp@sanger.ac.uk>

This file is part of cgpBattenberg.

cgpBattenberg is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

1. The usage of a range of years within a copyright statement contained within
this distribution should be interpreted as being equivalent to a list of years
including the first and last year specified and all consecutive years between
them. For example, a copyright statement that reads 'Copyright (c) 2005, 2007-
2009, 2011-2012' should be interpreted as being identical to a statement that
reads 'Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012' and a copyright
statement that reads "Copyright (c) 2005-2012' should be interpreted as being
identical to a statement that reads 'Copyright (c) 2005, 2006, 2007, 2008,
2009, 2010, 2011, 2012'."
```

<!-- Travis -->
[travis-base]: https://travis-ci.org/cancerit/cgpBattenberg
[travis-master]: https://travis-ci.org/cancerit/cgpBattenberg.svg?branch=master
[travis-develop]: https://travis-ci.org/cancerit/cgpBattenberg.svg?branch=dev

<!-- refs -->
[pcap-core-rel]: https://github.com/cancerit/PCAP-core/releases
[allele-count-rel]: https://github.com/cancerit/alleleCount/releases
[cgpvcf-rel]: https://github.com/cancerit/cgpVcf/releases
[impute-exe]: https://mathgen.stats.ox.ac.uk/impute/impute_v2.html
[bb-ref]: https://github.com/Wedge-Oxford/battenberg#required-reference-files
[cgpbb-wiki]: https://github.com/cancerit/cgpBattenberg/wiki
[bb-repo]: https://github.com/Wedge-Oxford/battenberg
[bb-ver-link]: https://github.com/Wedge-Oxford/battenberg/releases/tag/v2.2.5
[ds-cgpwgs-git]: https://github.com/cancerit/dockstore-cgpwgs

<!-- Quay.io -->
[quay-status]: https://quay.io/repository/wtsicgp/cgpbattenberg/status
[quay-repo]: https://quay.io/repository/wtsicgp/cgpbattenberg
[quay-builds]: https://quay.io/repository/wtsicgp/cgpbattenberg?tab=builds

