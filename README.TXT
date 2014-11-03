##########LICENCE##########
# Copyright (c) 2014 Genome Research Ltd.
#
# Author: Cancer Genome Project cgpit@sanger.ac.uk
#
# This file is part of battenberg.
#
# battenberg is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
##########LICENCE##########

Battenberg - An R program to detect subclonality and copy number in matched NGS data

####################
##	Installation	##
####################

	SEE INSTALL.TXT

################################
##	Program Run Instructions	##
################################

##Prerequisites-

	- Impute2 executable (https://mathgen.stats.ox.ac.uk/impute/impute_v2.html)
	******Any impute related data for download
	- BWA Mapped, indexed, duplicate marked/removed bam files, for both a matched normal and tumour sample
	- Reference.fasta and index
	- A file containing a list of contigs in the reference .fai to ignore
  - Directory containing the 1000 genomes allele and loci data (see File Formats below).
  - Impute info file impute_info.txt (see File Formats below)
  - Prob loci file probloci.txt (see File Formats below)


####################
##	File Formats	##
####################

##1000 Genomes data

****DW input required

##Impute info file

****DW input required

##Prob loci file

****DW input required


