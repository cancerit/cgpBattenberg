package Sanger::CGP::Battenberg::Implement;

##########LICENCE##########
# Copyright (c) 2014 Genome Research Ltd.
#
# Author: Cancer Genome Project cgpit@sanger.ac.uk
#
# This file is part of cgpBattenberg.
#
# cgpBattenberg is free software: you can redistribute it and/or modify it under
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


use strict;

use warnings FATAL => 'all';
use autodie qw(:all);
use Const::Fast qw(const);
use Cwd qw(abs_path getcwd);
use File::Basename;
use File::Spec;
use File::Which qw(which);
use File::Path qw(make_path remove_tree);
use File::Temp qw(tempfile);
use File::Copy qw(copy move);
use FindBin qw($Bin);
use List::Util qw(first);
use Vcf;
use POSIX qw(ceil);
use List::MoreUtils qw(uniq);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

use File::ShareDir qw(module_dir);

use Sanger::CGP::Battenberg;

use PCAP::Threaded;
use PCAP::Bam;

use Data::Dumper;

const my $ALLELE_COUNT_CMD => q{ -l %s -b %s -o %s -m %d -r %s};
const my $RUN_FUNC => q{ -e '%s %s'};
const my $SOURCE => q{source("%s"); };
const my $LIBRARY => q{library(Battenberg);};
const my $RUN_BAF_LOG => q{ getBAFsAndLogRs(tumourAlleleCountsFile.prefix="%s", normalAlleleCountsFile.prefix="%s", figuresFile.prefix="%s", BAFnormalFile="%s", BAFmutantFile="%s", logRnormalFile="%s", logRmutantFile="%s", combinedAlleleCountsFile="%s", %s, g1000file.prefix="%s", minCounts=%s, samplename="%s", seed=%s) };
const my $GC_CORRECT => q{ gc.correct.wgs(Tumour_LogR_file="%s", outfile="%s", correlations_outfile="%s", gc_content_file_prefix="%s", %s) };
const my $IMPUTE_FROM_AF => q{ generate.impute.input.wgs(chrom=%s, tumour.allele.counts.file="%s", normal.allele.counts.file="%s", output.file="%s", imputeinfofile="%s", is.male="%s", problemLociFile="%s", useLociFile=%s, heterozygousFilter=%s) };
const my $RUN_IMPUTE => q{ run.impute(inputfile="%s", outputfile.prefix="%s", is.male="%s", imputeinfofile="%s", impute.exe="%s", region.size=%s, chrom=%s, seed=%s)};
const my $COMBINE_IMPUTE => q{ combine.impute.output(inputfile.prefix="%s", outputfile="%s", is.male="%s", imputeinfofile="%s", region.size=%s, chrom=%s ) };
const my $HAPLOTYPE_BAF => q{ GetChromosomeBAFs(chrom="%s", SNP_file="%s", haplotypeFile="%s", samplename="%s", outfile="%s", %s, minCounts=%s) };
const my $PLOT_HAPLOTYPE_BAFS => q{ plot.haplotype.data(haplotyped.baf.file="%s", imageFileName="%s", samplename="%s", chrom=%s, %s) };
const my $COMBINE_BAFS => q{ combine.baf.files(inputfile.prefix="%s_chr", inputfile.postfix="_heterozygousMutBAFs_haplotyped.txt", outputfile="%s", no.chrs=%s) };
const my $SEGMENT_PHASED => q{ segment.baf.phased(samplename="%s", inputfile="%s", outputfile="%s", gamma=%s, phasegamma=%s, kmin=%s, phasekmin=%s, calc_seg_baf_option=%s ) };
const my $FIT_COPY_NUMBER => q{ fit.copy.number(samplename="%s", outputfile.prefix="%s_", inputfile.baf.segmented="%s", inputfile.baf="%s", inputfile.logr="%s", dist_choice=%s, ascat_dist_choice=%s, min.ploidy=%s, max.ploidy=%s, min.rho=%s, max.rho=%s, min.goodness=%s, uninformative_BAF_threshold=%s, gamma_param=%s %s) };
const my $CALL_SUBCLONES => q{ callSubclones(sample.name="%s", baf.segmented.file="%s", logr.file="%s", rho.psi.file="%s", output.file="%s", output.figures.prefix="%s", output.gw.figures.prefix="%s", masking_output_file="%s", sv_breakpoints_file="%s", %s, gamma=%s, segmentation.gamma=%s, siglevel=%s, maxdist=%s, noperms=%s, seed=%s, calc_seg_baf_option=%s) };
const my $GET_CHROM_NAMES => q{ get.chrom.names("%s", "%s")};
const my $ALLELE_COUNT_OUTPUT => q{%s_alleleFrequencies_chr%d.txt};
const my $ALLELE_COUNT_PREFIX => q{%s_alleleFrequencies_chr};
const my $ALLELE_LOCI_NAME => q{1000genomesloci2012_chr%d.txt};
const my $ALLELE_COUNT_TAR => q{%s_allelecounts.tar.gz};
const my $ALLELE_COUNT_DIR => q{%s_allelecounts};
const my $ALLELE_COUNT_SCRIPT => q{alleleCounter};
const my $NORMAL_CONTAMINATION_FILE => q{%s_normal_contamination.txt};
const my $RHO_PSI_FILE => q{%s_rho_and_psi.txt};
const my $SUBCLONE_PNG_OUTPUT => q{%s_subclones_chr%s.png};
const my $SUBCLONE_PREFIX => q{%s_subclones_chr};
const my $SUBCLONE_TAR => q{%s_subclones.tar.gz};
const my $SUBCLONE_DIR => q{%s_subclones};
const my $PROFILE_PNG => q{%s_BattenbergProfile};
const my $SEG_MASKING_TXT => q{%s_segment_masking_details.txt};
const my $HETDATA_OUTPUT => q{%s_chr%s_heterozygousData.png};
const my $HETDATA_TAR => q{%s_hetdata.tar.gz};
const my $HETDATA_DIR => q{%s_hetdata};
const my $HETBAFTXT_OUTPUT => q{%s_chr%s_heterozygousMutBAFs_haplotyped.txt};
const my $COMBINE_HETBAFTXT => q{%s_heterozygousMutBAFs_haplotyped.txt};
const my $HETBAFTXT_TAR => q{%s_hetbaf.tar.gz};
const my $HETBAFTXT_DIR => q{%s_hetbaf};
const my $OTHER_PNG_OUTPUT => q{%s_segment_chr%s.png};
const my $OTHER_PNG_TAR => q{%s_other.tar.gz};
const my $OTHER_PNG_DIR => q{%s_other};
const my $RAFSEG_PNG_OUTPUT => q{%s_RAFseg_chr%s.png};
const my $RAFSEG_PNG_TAR => q{%s_rafseg.tar.gz};
const my $RAFSEG_PNG_DIR => q{%s_rafseg};
const my $IMPUTE_INPUT_OUTPUT => q{%s_impute_input_chr%d.txt};
const my $IMPUTE_INPUT_TAR => q{%s_impute_input.tar.gz};
const my $IMPUTE_INPUT_DIR => q{%s_impute_input};
const my $IMPUTE_OUTPUT_OUTPUT => q{%s_impute_output_chr%d_allHaplotypeInfo.txt};
const my $IMPUTE_OUTPUT_TAR => q{%s_impute_output.tar.gz};
const my $IMPUTE_OUTPUT_DIR => q{%s_impute_output};
const my $SUNRISE_PNG	=> q{%s_distance.png};
const my $TUMOUR_PNG	=> q{%s.tumour.png};
const my $TUMOUR_CORRECTED_PNG	=> q{%s_Tumor.png};
const my $NORMAL_PNG	=> q{%s.germline.png};
const my $NORMAL_CORRECTED_PNG	=> q{%s_Germline.png};
const my $COPY_NO_PNG	=> q{%s_second_copynumberprofile.png};
const my $NON_ROUNDED_PNG	=> q{%s_second_nonroundedprofile.png};
const my $ALT_COPY_NO_ROUNDED_PNG => q{%s_copynumberprofile.png};
const my $ALT_NON_ROUNDED_PNG => q{%s_nonroundedprofile.png};
const my $SECOND_DISTANCE_PNG => q{%s_second_distance.png};
const my $PROFILE_AVERAGE_PNG => q{%s_BattenbergProfile_average.png};
const my $PROFILE_SUBCLONES_PNG => q{%s_BattenbergProfile_subclones.png};
const my $BAF_SEGMENT_TXT => q{%s.BAFsegmented.txt};
const my $LOGR_SEGMENT_TXT => q{%s.logRsegmented.txt};
const my $SUBCLONES_TXT => q{%s_subclones.txt};
const my $CELLULARITY_PLOIDY_TXT => q{%s_cellularity_ploidy.txt};
const my $SEGMENT_VCF => q{%s_logR_Baf_segmented.vcf};
const my $SEGMENT_VCF_GZ => q{%s_logR_Baf_segmented.vcf.gz};
const my $SEGMENT_VCF_TABIX => q{%s_logR_Baf_segmented.vcf.gz.tbi};
const my $CN_VCF => q{%s_battenberg_cn.vcf};
const my $CN_VCF_GZ => q{%s_battenberg_cn.vcf.gz};
const my $CN_VCF_TABIX => q{%s_battenberg_cn.vcf.gz.tbi};
const my $STATUS_TXT => q{%s_copynumber_solution_status.txt};
const my $SEED_TXT => q{%s_battenberg_seed.txt};
const my $GC_CORRECTED_TAB => q{%s_mutantLogR_gcCorrected.tab};
const my $GC_WINDOW_CORR_TXT => q{%s_GCwindowCorrelations.txt};
const my $GC_CORRECT_PREFIX => q{1000_genomes_GC_corr_chr_};
const my $NORMAL_BAF_TAB => q{%s_normalBAF.tab};
const my $TUMOUR_BAF_TAB => q{%s_mutantBAF.tab};
const my $NORMAL_LOGR_TAB => q{%s_normalLogR.tab};
const my $TUMOUR_LOGR_TAB => q{%s_mutantLogR.tab};
const my $ALLELECOUNTS_TAB => q{%s_alleleCounts.tab};

const my $RSCRIPT => q{Rscript};
const my $IMPUTE_EXE => q{impute2};

const my $IMPUTE_INPUT => q{%s_impute_input_chr%d_*K.*};
const my $IMPUTE_OUTPUT => q{%s_impute_output_chr%d.txt};

const my $ALLELE_COUNT_PARA => ' -b %s -o %s -l %s ';

#thousand genomes loci files
const my $ONEKGEN_ALLELE => q{1000genomesAlleles2012_chr};
const my $ONEKGEN_LOCI_FILE_PATTERN => q{1000genomesloci2012_chr*.txt};
const my $ONEKGEN_LOCI_FILE_REGEX => q{1000genomesloci2012_chr(\w+).txt};
const my $SPLIT_LOCI_GLOB => q{1000genomesloci2012_chr*_split%d.txt};
const my $SPLIT_LOCI_ALL_GLOB => q{1000genomesloci2012_chr*_split*.txt};
const my $SPLIT_LOCI_REGEX => q{1000genomesloci2012_chr(\w+)_split};
const my $SPLIT_ALLELE_COUNT_OUTPUT => q{%s_alleleFrequencies_chr%d_split%d.txt};
const my $SPLIT_ALLELE_COUNT_OUTPUT_GLOB => q{*_alleleFrequencies_chr*_split*.txt};
const my $SPLIT_ALLELE_COUNT_OUTPUT_REGEX => q{(.*)_alleleFrequencies_chr(\w+)_split(\d+).txt};

#runbaflog
const my $MIN_NORMAL_DEPTH => 10;

#imputefromaf
const my $HETEROZYGOUS_FILTER => 0.1;

#runimpute
const my $IMPUTE_REGION_SIZE => 5000000;

#segmentphased
const my $SEGMENTATION_GAMMA => 10;
const my $PHASING_GAMMA => 1;
const my $KMIN => 3;
const my $PHASE_KMIN => 3;

#callsubclones
const my $SIGLEVEL => 0.05;
const my $MAXDIST => 0.01;
const my $NOPERMS => 1000;

const my @BATTENBERG_RESULT_FILES => qw(
																					%s_Tumor.png
																					%s_Germline.png
																					%s_distance
																					%s_second_copynumberprofile
																					%s_second_nonroundedprofile
																					%s_nonroundedprofile
																					%s_second_distance
																					%s.BAFsegmented.txt
																					%s.logRsegmented.txt
																				);

sub get_mod_path {
  my $mod_path = dirname(abs_path($0)).'/../share';
  $mod_path = module_dir('Sanger::CGP::Battenberg::Implement') unless(-e File::Spec->catdir($mod_path, 'battenberg'));
  return $mod_path;
}

sub prepare {
  my $options = shift;
  if(exists $options->{'allele-counts'} && defined $options->{'allele-counts'}) {
    $options->{'tumour_name'} = $options->{'tumbam'};
    $options->{'normal_name'} = $options->{'normbam'};
  }
  else {
    $options->{'tumbam'} = File::Spec->rel2abs($options->{'tumbam'});
    $options->{'normbam'} = File::Spec->rel2abs($options->{'normbam'});
    $options->{'tumour_name'} = (PCAP::Bam::sample_name($options->{'tumbam'}, 1))[0];
    $options->{'normal_name'} = (PCAP::Bam::sample_name($options->{'normbam'}, 1))[0];
  }
	$options->{'mod_path'} = get_mod_path();
	$options->{'bat_path'} = File::Spec->catdir($options->{'mod_path'}, 'battenberg');
	$options->{'tmp'} = File::Spec->rel2abs($options->{'tmp'});
  return 1;
}

sub battenberg_splitlocifiles {
	my $options = shift;
	my $tmp = $options->{'tmp'};
	my $thou_gen_loc = File::Spec->rel2abs($options->{'1kgenloc'});
  my $requested_num_loci_files = $options->{'num_loci_files'};
  my $ignored_contigs = $options->{'ignored_contigs'};

	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  #No need to split if allele-counts parameter is set
  #(reading pre-calcuated allelecounter files)
	if(exists $options->{'allele-counts'} && defined $options->{'allele-counts'}) {
    return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
  }

  my $num_loci_files = $requested_num_loci_files;

  #Find the number of loci (lines) in each file and the grand total
  my ($loci_per_file, $total_loci_remaining) = _find_num_loci_per_file($thou_gen_loc, $ignored_contigs);

  #Average number of loci (lines) in a split file
  my $num_loci_per_split = $total_loci_remaining / $num_loci_files;

  my $total_number_of_split_files_created = 0;

  foreach my $file (keys %$loci_per_file) {

    #Check if we are at the last file
    my $number_of_split_files_for_chr;
    if ($total_loci_remaining == $loci_per_file->{$file}) {
      #Set last file to be what is left
      $number_of_split_files_for_chr = $requested_num_loci_files - $total_number_of_split_files_created;
    } else {
      $number_of_split_files_for_chr = int($loci_per_file->{$file} / $num_loci_per_split);
    }

    #Must have each file at least once
    if ($number_of_split_files_for_chr == 0) {
      $number_of_split_files_for_chr = 1;
    }
    #print "file=$file " . $loci_per_file->{$file} . " $number_of_split_files_for_chr $total_number_of_split_files_created\n";

    #Write split loci files to the tmpdir
    _create_split_files($file, $loci_per_file->{$file}, $number_of_split_files_for_chr, $tmp, $total_number_of_split_files_created);

    #Keep track of the total number of files created
    $total_number_of_split_files_created += $number_of_split_files_for_chr;

    #Need to keep refining the answer to avoid problems with rounding
    if ($num_loci_files > $number_of_split_files_for_chr) {
      $num_loci_per_split = ($total_loci_remaining -= $loci_per_file->{$file}) / ($num_loci_files -= $number_of_split_files_for_chr);
    }
  }

  if ($total_number_of_split_files_created != $requested_num_loci_files) {
    die("The number of loci files created $total_number_of_split_files_created is not the same as the number requested $num_loci_files\n");
  }
  return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub battenberg_allelecount {
	# uncoverable subroutine
	my ($index_in, $options) = @_;
	return 1 if(exists $options->{'index'} && $index_in != $options->{'index'});
	my $tmp = $options->{'tmp'};

  #Unpack the allele counts tar file if it exists and return
	if(exists $options->{'allele-counts'} && defined $options->{'allele-counts'}) {
	  unless ($index_in == 1) {
	    return 1;
	  }
	  # expand the file and put the data in the expected locations
	  my $command = sprintf 'tar -C %s -zxf %s', $tmp, $options->{'allele-counts'};

	  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index_in);

	  # then they have to be renamed:
	  # %s_alleleFrequencies_chr%d.txt
	  my $sampdir = File::Spec->catdir($tmp, $options->{'tumour_name'});
	  opendir(my $dh, $sampdir) || die "ERROR opening $sampdir: $!";
	  while(my $thing = readdir $dh) {
  	  if($thing =~ m/^(.*)[.]([^.]+)[.]tsv$/) {
  	    my $samp = $1;
  	    my $chr = $2;
  	    my $moved = File::Spec->catfile($tmp, sprintf $ALLELE_COUNT_OUTPUT, $samp, $chr);
  	    my $orig = File::Spec->catfile($sampdir, $thing);
  	    move($orig, $moved) || die "Unable to move file $orig -> $moved";
  	  }
  	}
	  remove_tree(File::Spec->catdir($tmp, $options->{'tumour_name'}));
    return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index_in);
  }

  #####
  #Run allelecounter program
  #####
  my $reference = $options->{'reference'};
  my $mbq =  $options->{'mbq'};
  my $num_loci_files = $options->{'num_loci_files'};

  my @indices = limited_indices($options, $index_in, $num_loci_files*2);
  for my $index(@indices) {
    my $sample_name = $options->{'tumour_name'};
    my $sample_bam = $options->{'tumbam'};
    my $file_index = $index;
    if($index > $num_loci_files) {
      $sample_name = $options->{'normal_name'};
      $sample_bam = $options->{'normbam'};
      $file_index = $index - $num_loci_files;
    }

    #Skip if this job has already been successfully run
    next if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);

    my $loci_file;
    my $alleleCountOut;
    my $split_name = File::Spec->catfile($tmp, sprintf $SPLIT_LOCI_GLOB, $file_index);
    my @split_files = glob($split_name);

    #Found split file
    if (@split_files > 0) {
      #Should only find one
      if (@split_files > 1) {
        die("Found more than one file for index $index\n");
      }
      $loci_file = shift @split_files;
      my ($chr) = $loci_file =~ /$SPLIT_LOCI_REGEX/;
      $alleleCountOut = File::Spec->catfile($tmp, sprintf $SPLIT_ALLELE_COUNT_OUTPUT, $sample_name, $chr, $file_index);
    } else {
      my $k_gen_loc = File::Spec->rel2abs($options->{'1kgenloc'});

      #Get the appropriate loci file based on the index into the contigs array rather than the number
      #of requested contigs
      my $loci_names_to_index = _lociNameMap($k_gen_loc);
      my $contigs = read_contigs_from_file_with_ignore($options->{'reference'},$options->{'ignored_contigs'});
      my $this_contig = $contigs->[$file_index-1];
      my $file_index = $loci_names_to_index->{$this_contig};

      $loci_file = File::Spec->catfile($k_gen_loc,sprintf($ALLELE_LOCI_NAME,$file_index));
      $alleleCountOut = File::Spec->rel2abs(File::Spec->catfile($tmp,sprintf($ALLELE_COUNT_OUTPUT,$sample_name,$file_index)));
    }

    PCAP::Cli::file_for_reading('1k-genome-loci-file',$loci_file);
    my $command = _which($ALLELE_COUNT_SCRIPT) || die "Unable to find $ALLELE_COUNT_SCRIPT in path";

    $command .= sprintf($ALLELE_COUNT_CMD,
                        $loci_file,
                        $sample_bam,
                        $alleleCountOut,
                        $mbq,
                        $reference,
                       );
    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
  }
}

sub battenberg_concat_allelecount {
  my ($options) = @_;
  my $tmp = $options->{'tmp'};

  #Skip if already done
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  #No need to concatenate if allele-counts parameter is set
  #(reading pre-calcuated allelecounter files)
	if(exists $options->{'allele-counts'} && defined $options->{'allele-counts'}) {
    return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
  }

  #Find all the loci split files
  my $split_name = File::Spec->catfile($tmp, sprintf $SPLIT_LOCI_ALL_GLOB,);
  my @split_files = glob($split_name);

  #Return if there are no split files, so nothing to concatentate
  unless (@split_files) {
    return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
  }

  #Find how many split files there should be for each chromosome
  my %loci_files;
  foreach my $loci_file (@split_files) {
    my ($chr) = $loci_file =~ /$SPLIT_LOCI_REGEX/;
    $loci_files{$chr}++;
  }

  #Find all the allele count files
  my $allele_count_filenames = File::Spec->catfile($tmp, $SPLIT_ALLELE_COUNT_OUTPUT_GLOB);
  my @allele_count_files = glob($allele_count_filenames);

  my $sample_jobs;
  foreach my $split_file (@allele_count_files) {
    my ($basename, $path) = fileparse($split_file);

    #Grab the sample name, chromosome and order from the filename
    my ($sname, $chr, $file_cnt) = $basename =~ /$SPLIT_ALLELE_COUNT_OUTPUT_REGEX/;

    #All the files with the same sample name and chromosome need concatenating together
    my $split_file_cnt;
    %$split_file_cnt = ('order' => $file_cnt,
                        'file' => $split_file);
    push @{$sample_jobs->{$sname}{$chr}}, $split_file_cnt;
  }

  foreach my $sname (keys %$sample_jobs) {
    foreach my $contig_cnt (keys %{$sample_jobs->{$sname}}) {

      #Check we have the correct number of files
      if ($loci_files{$contig_cnt} != @{$sample_jobs->{$sname}{$contig_cnt}}) {
        die('ERROR: The number of allele count files (' . @{$sample_jobs->{$sname}{$contig_cnt}} . ') does not equal the number of loci files (' . $loci_files{$contig_cnt} . ')');
      }

      #Create allele count file for each sample and chromosome
      my $alleleCountOut = File::Spec->rel2abs(File::Spec->catfile($tmp,sprintf($ALLELE_COUNT_OUTPUT,$sname,$contig_cnt)));
      #Make sure any existing file is deleted before appending
      unlink $alleleCountOut if (-e $alleleCountOut);

      open my $out_fh, ">> $alleleCountOut" or die "Unable to open $alleleCountOut for writing";

      my $first_file = 1;
      foreach my $split_file (sort {$a->{'order'} <=> $b->{'order'}} @{$sample_jobs->{$sname}{$contig_cnt}}) {
        my $filename =  $split_file->{'file'};
        open my $in_fh, '<', $filename or die "Unable to open $filename";
        while (my $line = <$in_fh>) {
          #print comments of the first file, then skip for all other files
          if ($line =~ /^#/) {
            if ($first_file) {
              print $out_fh $line;
            }
          } else {
            #print rest of file
            print $out_fh $line;
          }
        }
        close $in_fh;
        $first_file = 0;
      }
      close $out_fh;
    }
  }
  return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub battenberg_runbaflog{
	# uncoverable subroutine
	my $options = shift;
	my $tmp = $options->{'tmp'};
	my $mod_path = $options->{'mod_path'};
	my $impute_info = $options->{'impute_info'};
	my $thou_gen_loc = File::Spec->rel2abs($options->{'1kgenloc'});
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

	#Use the share directory so we don't have to symlink anything.
	my $command = "cd $tmp; ";
	$command .= _which($RSCRIPT) || die "Unable to find $RSCRIPT in path";

  my $tumour_name = $options->{'tumour_name'};
  my $normal_name = $options->{'normal_name'};
  my $tumourAlleleCountsFile_prefix = sprintf $ALLELE_COUNT_PREFIX, $tumour_name;
  my $normalAlleleCountsFile_prefix = sprintf $ALLELE_COUNT_PREFIX, $normal_name;

  my $figuresFile_prefix = $tumour_name . '_';
  my $BAFnormalFile = sprintf $NORMAL_BAF_TAB, $tumour_name;
  my $BAFmutantFile = sprintf $TUMOUR_BAF_TAB, $tumour_name;
  my $logRnormalFile = sprintf $NORMAL_LOGR_TAB, $tumour_name;
  my $logRmutantFile = sprintf $TUMOUR_LOGR_TAB, $tumour_name;
  my $combinedAlleleCountsFile = sprintf $ALLELECOUNTS_TAB, $tumour_name;
  my $g1000file_prefix = File::Spec->catfile($thou_gen_loc, $ONEKGEN_ALLELE);
  my $minCounts = $MIN_NORMAL_DEPTH,
  my $samplename = $options->{'tumour_name'};
  my $seed = $options->{'seed'};

  my $chr_names = _get_chroms_as_string($impute_info, $options->{'is_male'});

#  my $source1 = sprintf $SOURCE, File::Spec->catfile($options->{'bat_path'}, 'prepare_wgs.R');
#  my $source2 = sprintf $SOURCE, File::Spec->catfile($options->{'bat_path'}, 'util.R');
#  my $source3 = sprintf $SOURCE, File::Spec->catfile($options->{'bat_path'}, 'clonal_ascat.R');

  my $function = sprintf $RUN_BAF_LOG,
                         $tumourAlleleCountsFile_prefix,
                         $normalAlleleCountsFile_prefix,
                         $figuresFile_prefix,
                         $BAFnormalFile,
                         $BAFmutantFile,
                         $logRnormalFile,
                         $logRmutantFile,
                         $combinedAlleleCountsFile,
                         "chr_names=as.vector(c($chr_names))",
                         $g1000file_prefix,
                         $minCounts,
                         $samplename,
                         $seed
                         ;
  #$command .= sprintf $RUN_FUNC, "$source1 $source2 $source3", $function;
  $command .= sprintf $RUN_FUNC, $LIBRARY, $function;

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);

	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub battenberg_gc_correct {
  # uncoverable subroutine
	my ($options) = @_;
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
	my $mod_path = $options->{'mod_path'};
	my $impute_info = $options->{'impute_info'};
	my $prob_loci = $options->{'prob_loci'};
  my $is_male = $options->{'is_male'};
  my $tumour_name = $options->{'tumour_name'};
  my $normal_name = $options->{'normal_name'};
  my $logRmutantFile = sprintf $TUMOUR_LOGR_TAB, $tumour_name;
  my $outfile = sprintf $GC_CORRECTED_TAB, $tumour_name;
  my $correlations_outfile = sprintf $GC_WINDOW_CORR_TXT, $tumour_name;
  my $gc_content_file_prefix = File::Spec->catfile($options->{'gc_correct_loc'}, $GC_CORRECT_PREFIX);
  my $chr_names = _get_chroms_as_string($impute_info, $options->{'is_male'});

	#Use the share directory so we don't have to symlink anything.
	my $command = "cd $tmp; ";
	$command .= _which($RSCRIPT) || die "Unable to find $RSCRIPT in path";

  #my $source1 = sprintf $SOURCE, File::Spec->catfile($options->{'bat_path'}, 'prepare_wgs.R');
  my $function = sprintf $GC_CORRECT,
    $logRmutantFile,
    $outfile,
    $correlations_outfile,
    $gc_content_file_prefix,
    "chrom_names=as.vector(c($chr_names))";

#  $command .= sprintf $RUN_FUNC, "$source1", $function;
  $command .= sprintf $RUN_FUNC, $LIBRARY, $function;

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub battenberg_imputefromaf{
	# uncoverable subroutine
	my ($index, $options) = @_;
	return 1 if(exists $options->{'index'} && $index != $options->{'index'});
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);
	my $mod_path = $options->{'mod_path'};
	my $impute_info = $options->{'impute_info'};
	my $prob_loci = $options->{'prob_loci'};
  my $is_male = $options->{'is_male'};
  my $tumour_name = $options->{'tumour_name'};
  my $normal_name = $options->{'normal_name'};

  my $tumour_allele_counts_file =  File::Spec->catfile(sprintf($ALLELE_COUNT_OUTPUT,$tumour_name, $index));
  my $normal_allele_counts_file =  File::Spec->catfile(sprintf($ALLELE_COUNT_OUTPUT,$normal_name, $index));

  my $output_file = sprintf $IMPUTE_INPUT_OUTPUT, $tumour_name, $index;
  my $use_loci_file = 'NA';
  my $heterozygousFilter = $HETEROZYGOUS_FILTER;

	#Use the share directory so we don't have to symlink anything.
	my $command = "cd $tmp; ";
	$command .= _which($RSCRIPT) || die "Unable to find $RSCRIPT in path";

#  my $source1 = sprintf $SOURCE, File::Spec->catfile($options->{'bat_path'}, 'prepare_wgs.R');
#  my $source2 = sprintf $SOURCE, File::Spec->catfile($options->{'bat_path'}, 'impute.R');
  my $function =sprintf $IMPUTE_FROM_AF,
    $index,
    $tumour_allele_counts_file,
    $normal_allele_counts_file,
    $output_file,
    $impute_info,
    $is_male,
    $prob_loci,
    $use_loci_file,
    $heterozygousFilter;

#  $command .= sprintf $RUN_FUNC, "$source1 $source2", $function;
  $command .= sprintf $RUN_FUNC, $LIBRARY, $function;

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
}

sub battenberg_runimpute{
	# uncoverable subroutine
	my ($index, $options) = @_;
	return 1 if(exists $options->{'index'} && $index != $options->{'index'});
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);
	my $mod_path = $options->{'mod_path'};
	my $impute_info = $options->{'impute_info'};
	my $impute_exe = _which($IMPUTE_EXE) || die "Unable to find $IMPUTE_EXE in path";

  my $tumour_name = $options->{'tumour_name'};
  my $normal_name = $options->{'normal_name'};
  my $input_file = sprintf $IMPUTE_INPUT_OUTPUT, $tumour_name, $index;
  my $output_file = sprintf $IMPUTE_OUTPUT, $tumour_name, $index;
  my $is_male = $options->{'is_male'};
  my $seed = $options->{'seed'};
  my $region_size = $IMPUTE_REGION_SIZE;

	#Use the share directory so we don't have to symlink anything.
	my $command = "cd $tmp; ";

	$command .= _which($RSCRIPT) || die "Unable to find $RSCRIPT in path";
#  my $source1 = sprintf $SOURCE, File::Spec->catfile($options->{'bat_path'}, 'impute.R');
  my $function =sprintf $RUN_IMPUTE,
    $input_file,
    $output_file,
    $is_male,
    $impute_info,
    $impute_exe,
    $region_size,
    $index,
    $seed;

 # $command .= sprintf $RUN_FUNC, "$source1", $function;
  $command .= sprintf $RUN_FUNC, $LIBRARY, $function;

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
}

sub battenberg_combineimpute{
	# uncoverable subroutine
	my ($index, $options) = @_;
	return 1 if(exists $options->{'index'} && $index != $options->{'index'});
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);
	my $mod_path = $options->{'mod_path'};
	my $impute_info = $options->{'impute_info'};

  my $tumour_name = $options->{'tumour_name'};
  my $normal_name = $options->{'normal_name'};
  my $input_file = sprintf $IMPUTE_OUTPUT, $tumour_name, $index;
  my $output_file = sprintf $IMPUTE_OUTPUT_OUTPUT, $tumour_name, $index;

  my $is_male = $options->{'is_male'};
  my $region_size = $IMPUTE_REGION_SIZE;

	#Use the share directory so we don't have to symlink anything.
	my $command = "cd $tmp; ";
	$command .= _which($RSCRIPT) || die "Unable to find $RSCRIPT in path";

#  my $source1 = sprintf $SOURCE, File::Spec->catfile($options->{'bat_path'}, 'impute.R');
#  my $source2 = sprintf $SOURCE, File::Spec->catfile($options->{'bat_path'}, 'util.R');
  my $function =sprintf $COMBINE_IMPUTE,
    $input_file,
    $output_file,
    $is_male,
    $impute_info,
    $region_size,
    $index;

#  $command .= sprintf $RUN_FUNC, "$source1 $source2", $function;
  $command .= sprintf $RUN_FUNC, $LIBRARY, $function;

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
}

sub battenberg_haplotypebaf{
	# uncoverable subroutine
	my ($index, $options) = @_;
	return 1 if(exists $options->{'index'} && $index != $options->{'index'});
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);
	my $mod_path = $options->{'mod_path'};
	my $impute_info = $options->{'impute_info'};

  my $tumour_name = $options->{'tumour_name'};
  my $normal_name = $options->{'normal_name'};
  my $snp_file =  File::Spec->catfile(sprintf($ALLELE_COUNT_OUTPUT,$tumour_name, $index));

  my $haplotype_file = sprintf $IMPUTE_OUTPUT_OUTPUT, $tumour_name, $index;
  my $outfile = sprintf $HETBAFTXT_OUTPUT, $tumour_name, $index;
  my $min_counts = $MIN_NORMAL_DEPTH;
  my $chr_names = _get_chroms_as_string($impute_info, $options->{'is_male'});

	my $command = "cd $tmp; ";
	$command .= _which($RSCRIPT) || die "Unable to find $RSCRIPT in path";

#  my $source1 = sprintf $SOURCE, File::Spec->catfile($options->{'bat_path'}, 'haplotype.R');
  my $function =sprintf $HAPLOTYPE_BAF,
    $index,
    $snp_file,
    $haplotype_file,
    $tumour_name, $outfile,
    "chr_names=as.vector(c($chr_names))",
     $min_counts;

#  $command .= sprintf $RUN_FUNC, "$source1", $function;
  $command .= sprintf $RUN_FUNC, $LIBRARY, $function;

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
}

sub battenberg_postbafcleanup{
	# uncoverable subroutine
	my ($index, $options) = @_;
	return 1 if(exists $options->{'index'} && $index != $options->{'index'});
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);
	my $mod_path = $options->{'mod_path'};

	my $impute_in = sprintf($IMPUTE_INPUT,$options->{'tumour_name'},$index);
	my $impute_out = sprintf($IMPUTE_OUTPUT,$options->{'tumour_name'},$index);

	my $command = "rm -f $impute_in $impute_out";

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
}

sub battenberg_plothaplotypes{
	# uncoverable subroutine
	my ($index, $options) = @_;
	return 1 if(exists $options->{'index'} && $index != $options->{'index'});
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);
	my $mod_path = $options->{'mod_path'};
	my $impute_info = $options->{'impute_info'};
  my $tumour_name = $options->{'tumour_name'};
  my $normal_name = $options->{'normal_name'};
  my $haplotype_file = sprintf $HETBAFTXT_OUTPUT, $tumour_name, $index;
  my $image_file = sprintf $HETDATA_OUTPUT, $tumour_name, $index;
  my $chr_names = _get_chroms_as_string($impute_info, $options->{'is_male'});

	my $command = "cd $tmp; ";
	$command .= _which($RSCRIPT) || die "Unable to find $RSCRIPT in path";

#  my $source1 = sprintf $SOURCE, File::Spec->catfile($options->{'bat_path'}, 'haplotype.R');
#  my $source2 = sprintf $SOURCE, File::Spec->catfile($options->{'bat_path'}, 'plotting.R');
  my $function =sprintf $PLOT_HAPLOTYPE_BAFS,
    $haplotype_file,
    $image_file,
    $tumour_name,
    $index,
    "chr_names=as.vector(c($chr_names))";

#  $command .= sprintf $RUN_FUNC, "$source1 $source2", $function;
  $command .= sprintf $RUN_FUNC, $LIBRARY, $function;

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
}

sub battenberg_combinebafs{
	# uncoverable subroutine
	my $options = shift;
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
	my $mod_path = $options->{'mod_path'};
	my $impute_info = $options->{'impute_info'};

  my $tumour_name = $options->{'tumour_name'};
  my $normal_name = $options->{'normal_name'};
  my $output_file = sprintf $COMBINE_HETBAFTXT, $tumour_name;
  my $chroms = _get_chrom_names($impute_info, $options->{'is_male'});
  my $num_chrs = @$chroms;

	my $command = "cd $tmp; ";
	$command .= _which($RSCRIPT) || die "Unable to find $RSCRIPT in path";

  #my $source1 = sprintf $SOURCE, File::Spec->catfile($options->{'bat_path'}, 'haplotype.R');
  #my $source2 = sprintf $SOURCE, File::Spec->catfile($options->{'bat_path'}, 'util.R');
  my $function =sprintf $COMBINE_BAFS,
    $tumour_name,
    $output_file,
    $num_chrs;

#  $command .= sprintf $RUN_FUNC, "$source1 $source2", $function;
  $command .= sprintf $RUN_FUNC, $LIBRARY, $function;

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub battenberg_segmentphased{
	# uncoverable subroutine
	my $options = shift;
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  my $tumour_name = $options->{'tumour_name'};
  my $normal_name = $options->{'normal_name'};
  my $input_file = sprintf $COMBINE_HETBAFTXT, $tumour_name;
  my $output_file = sprintf $BAF_SEGMENT_TXT, $tumour_name;
  my $gamma = $SEGMENTATION_GAMMA;
  my $phasegamma = $PHASING_GAMMA;
  my $kmin = $KMIN;
  my $phasekmin = $PHASE_KMIN;
  my $calc_seg_baf_option = $options->{'calc_seg_baf_option'};

	my $command = "cd $tmp; ";
	$command .= _which($RSCRIPT) || die "Unable to find $RSCRIPT in path";

 # my $source1 = sprintf $SOURCE, File::Spec->catfile($options->{'bat_path'}, 'segmentation.R');
 # my $source2 = sprintf $SOURCE, File::Spec->catfile($options->{'bat_path'}, 'fastPCF.R');
 # my $source3 = sprintf $SOURCE, File::Spec->catfile($options->{'bat_path'}, 'plotting.R');

  my $function = sprintf $SEGMENT_PHASED,
    $tumour_name,
    $input_file,
    $output_file,
    $gamma,
    $phasegamma,
    $kmin,
    $phasekmin,
    $calc_seg_baf_option;

#  $command .= sprintf $RUN_FUNC, "$source1 $source2 $source3", $function;
  $command .= sprintf $RUN_FUNC, $LIBRARY, $function;

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);

}

sub battenberg_fitcopyno{
	# uncoverable subroutine
	my $options = shift;
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
  my $tumour_name = $options->{'tumour_name'};
  my $normal_name = $options->{'normal_name'};
  my $BAFsegmentedFile = sprintf $BAF_SEGMENT_TXT, $tumour_name;
  my $BAFmutantFile = sprintf $TUMOUR_BAF_TAB, $tumour_name;
  my $input_file_logr = sprintf $TUMOUR_LOGR_TAB, $tumour_name ;

	my $command = "cd $tmp; ";
	$command .= _which($RSCRIPT) || die "Unable to find $RSCRIPT in path";

#  my $source1 = sprintf $SOURCE, File::Spec->catfile($options->{'bat_path'}, 'fitcopynumber.R');
#  my $source2 = sprintf $SOURCE, File::Spec->catfile($options->{'bat_path'}, 'util.R');
#  my $source3 = sprintf $SOURCE, File::Spec->catfile($options->{'bat_path'}, 'clonal_ascat.R');
#  my $source4 = sprintf $SOURCE, File::Spec->catfile($options->{'bat_path'}, 'orderEdges.R');
#  my $source5 = sprintf $SOURCE, File::Spec->catfile($options->{'bat_path'}, 'plotting.R');

 if(exists($options->{'new_chr'}) && defined($options->{'new_chr'})){
    my ($rho, $psi) = _calc_rho_psi_from_subclones_file($options);
    #put rho and psi in the params.
    $options->{'psi'} = $psi;
    $options->{'rho'} = $rho;
  }
  my $rho_psi_cmds = q{};
  if (exists ($options->{'rho'}) && defined($options->{'rho'})){
    #Sanity checks for rho only. Truncate to 5 decimal places for max comparison to allow for rounding issues
    my $round_rho = sprintf("%.5f", $options->{'rho'});
    if ($options->{'rho'} < $options->{'min_rho'} || $round_rho > $options->{'max_rho'}) {
      die("Invalid value for rho (" . $round_rho . "). This should be in the range specified by min_rho (" . $options->{'min_rho'} . ") and max_rho (" . $options->{'max_rho'} . ")");
    }
    $rho_psi_cmds .= ', use_preset_rho_psi=T';
    $rho_psi_cmds .= ', preset_rho='. $options->{'rho'};
    $rho_psi_cmds .= ', preset_psi=' . $options->{'psi'};
  }
  my $function = sprintf $FIT_COPY_NUMBER,
    $tumour_name,
    $tumour_name,
    $BAFsegmentedFile,
    $BAFmutantFile,
    $input_file_logr,
    $options->{'clonality_dist'},
    $options->{'ascat_dist'},
    $options->{'min_ploidy'},
    $options->{'max_ploidy'},
    $options->{'min_rho'},
    $options->{'max_rho'},
    $options->{'min_goodness'},
    $options->{'balanced_thresh'},
    $options->{'plat_gamma'},
    $rho_psi_cmds;

  $command .=  sprintf $RUN_FUNC, $LIBRARY, $function;

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub _calc_rho_psi_from_subclones_file {
  my ($options) = @_;

  my $tmp = $options->{'tmp'};
  my $new_chr = $options->{'new_chr'};
  my $new_pos = $options->{'new_pos'};
  my $min_cn = $options->{'new_min_cn'};
  my $maj_cn= $options->{'new_maj_cn'};

  #Get the subclones file.
  my $subcl_txt = File::Spec->catfile($tmp,sprintf("$SUBCLONES_TXT.gz",$options->{'tumour_name'}));
  die("Could not locate subclones.txt file '$subcl_txt'.") unless(-e $subcl_txt);
  #open and find the appropriate line
  my $ref_baf;
  my $log_r_ref;
  my $FH = IO::Uncompress::Gunzip->new($subcl_txt) or die "IO::Uncompress::Gunzip failed: $GunzipError\n";

  while(<$FH>){
    my $line = $_;
    next if($line =~ m/\s*chr/);
    chomp($line);
    my ($chr,$start,$stop,$BAF,$unused,$LogR,undef) = split(/\t/,$line);
    if($chr eq $new_chr && $new_pos >= $start && $new_pos <= $stop) {
      $ref_baf = $BAF;
      $log_r_ref = $LogR;
      warn Dumper($chr,$start,$stop,$ref_baf,$log_r_ref);
      last;
    }

  }
  close($FH) || die("Error trying to close subclones file '$subcl_txt' :$!");

  die("Unable to find location ($new_chr:$new_pos) in subclones file $subcl_txt") unless (defined $ref_baf);
  #refChromosome and refPosition should be used to obtain refBAF for the segment containing this chr & position,
  #i.e. from the subclones file find the row that has refChromosome in column 1,
  #refPosition between col2 and col3 and take the value in col 4

  #use the information to calculate psi and rho.
  #The conversion formulae you need are:
  #rho = (2*refBAF-1)/(2*refBAF-refBAF*(refMajor+refMinor)-1+refMajor)
  #psi = (rho*(refMajor+refMinor)+2-2*rho)/(2^(LogRref/gamma))
  #psit = (psi-2*(1-rho))/rho
  my $gamma = 1;
  my $rho = (2 * $ref_baf - 1) / (2 * $ref_baf - $ref_baf * ($maj_cn+$min_cn) - 1 + $maj_cn);
  my $psi_tmp = ($rho * ($maj_cn + $min_cn) + 2 - 2*$rho)/(2**($log_r_ref/$gamma));
  my $psi = ($psi_tmp-2*(1-$rho))/$rho;

  #Remove the unrequired values from the params hash.
  delete $options->{'new_chr'};
  delete $options->{'new_pos'};
  delete $options->{'new_min_cn'};
  delete $options->{'new_maj_cn'};

  return($rho, $psi);
}

sub battenberg_callsubclones{
	# uncoverable subroutine
	my $options = shift;
	my $tmp = $options->{'tmp'};
	my $impute_info = $options->{'impute_info'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
  my $tumour_name = $options->{'tumour_name'};
  my $normal_name = $options->{'normal_name'};
  my $BAFsegmentedFile = sprintf $BAF_SEGMENT_TXT, $tumour_name;
  my $logr_file = sprintf $GC_CORRECTED_TAB, $tumour_name;
  my $rho_psi_file = sprintf $RHO_PSI_FILE,  $tumour_name;
  my $output_file = sprintf $SUBCLONES_TXT, $tumour_name;
  my $output_figures_prefix = sprintf $SUBCLONE_PREFIX, $tumour_name;
  my $output_gw_figures_prefix = sprintf $PROFILE_PNG, $tumour_name;
  my $masking_output_file = sprintf $SEG_MASKING_TXT, $tumour_name;
  my $sv_breakpoints_file="NA";
  my $chrom_names = _get_chroms_as_string($impute_info, $options->{'is_male'});
  my $segmentation_gamma="NA";
  my $siglevel = $SIGLEVEL;
  my $maxdist = $MAXDIST;
  my $noperms = $NOPERMS;
  my $calc_seg_baf_option = $options->{'calc_seg_baf_option'};
  my $seed = $options->{'seed'};

	my $command = "cd $tmp; ";
	$command .= _which($RSCRIPT) || die "Unable to find $RSCRIPT in path";

#  my $source1 = sprintf $SOURCE, File::Spec->catfile($options->{'bat_path'}, 'fitcopynumber.R');
#  my $source2 = sprintf $SOURCE, File::Spec->catfile($options->{'bat_path'}, 'util.R');
#  my $source3 = sprintf $SOURCE, File::Spec->catfile($options->{'bat_path'}, 'orderEdges.R');

  my $function = sprintf $CALL_SUBCLONES,
    $tumour_name,
    $BAFsegmentedFile,
    $logr_file,
    $rho_psi_file,
    $output_file,
    $output_figures_prefix,
    $output_gw_figures_prefix,
    $masking_output_file,
    $sv_breakpoints_file,
    "chr_names=as.vector(c($chrom_names))",
    $options->{'plat_gamma'},
    $segmentation_gamma,
    $siglevel,
    $maxdist,
    $noperms,
    $seed,
    $calc_seg_baf_option;

  $command .= sprintf $RUN_FUNC, $LIBRARY, $function;

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub battenberg_finalise{
	# uncoverable subroutine
	my $options = shift;
  my $tmp = $options->{'tmp'};

  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

 	my $outdir = $options->{'outdir'};
  my $tumour_name = $options->{'tumour_name'};
	my ($rho,$psi) = _extractRhoPsiFromFile($options);

	#Calculate normal contamination
	my $normal_contamination = 2*(1-$rho) /
												(2*(1-$rho) + $psi * $rho);
	$options->{'normc'} = $normal_contamination;
	_writeNormalContaminationToFile($options);

  #Write seed to file
  _writeSeedToFile($options);

	#Tarball allelcounts and copy to results folder
	if (PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), @{['allele_tar_gz',0]}) == 0){
		_zip_and_tar_fileset($options,
														$ALLELE_COUNT_OUTPUT,
														sprintf($ALLELE_COUNT_TAR,$options->{'normal_name'}),
														sprintf($ALLELE_COUNT_DIR,$options->{'normal_name'}),
														$tmp,
														$outdir,
														undef,
														$options->{'normal_name'});

		_zip_and_tar_fileset($options,
														$ALLELE_COUNT_OUTPUT,
														sprintf($ALLELE_COUNT_TAR,$options->{'tumour_name'}),
														sprintf($ALLELE_COUNT_DIR,$options->{'tumour_name'}),
														$tmp,
														$outdir,
														undef,
														$options->{'tumour_name'});

		PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), @{['allele_tar_gz',0]});
	}

	my $contigs = read_contigs_from_file_with_ignore($options->{'reference'},$options->{'ignored_contigs'});
	#Tarball subclones pngs and move to the results folder
	if(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), @{['subclone.tar.gz',0]}) == 0){
		_zip_and_tar_fileset($options,
														$SUBCLONE_PNG_OUTPUT,
														sprintf($SUBCLONE_TAR,$options->{'tumour_name'}),
														sprintf($SUBCLONE_DIR,$options->{'tumour_name'}),
														$tmp,
														$outdir,
														$contigs,
														$options->{'tumour_name'});
		PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), @{['subclone.tar.gz',0]});
	}

	#Tarball of by chromosome heterozygousData.png images
	if(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), @{['hetDataPngs.tar.gz',0]}) == 0){
		_zip_and_tar_fileset($options,
														$HETDATA_OUTPUT,
														sprintf($HETDATA_TAR,$options->{'tumour_name'}),
														sprintf($HETDATA_DIR,$options->{'tumour_name'}),
														$tmp,
														$outdir,
														$contigs,
														$options->{'tumour_name'});
		PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), @{['hetDataPngs.tar.gz',0]});
	}

	#het_mut_baf_txt tar.gz
	if(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), @{['hetBAFTxt.tar.gz',0]}) == 0){
		_zip_and_tar_fileset($options,
														$HETBAFTXT_OUTPUT,
														sprintf($HETBAFTXT_TAR,$options->{'tumour_name'}),
														sprintf($HETBAFTXT_DIR,$options->{'tumour_name'}),
														$tmp,
														$outdir,
														undef,
														$options->{'tumour_name'});
		PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), @{['hetBAFTxt.tar.gz',0]});
	}

	#Other pngs
	if(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), @{['other.tar.gz',0]}) == 0){
		_zip_and_tar_fileset($options,
														$OTHER_PNG_OUTPUT,
														sprintf($OTHER_PNG_TAR,$options->{'tumour_name'}),
														sprintf($OTHER_PNG_DIR,$options->{'tumour_name'}),
														$tmp,
														$outdir,
														$contigs,
														$options->{'tumour_name'});
		PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), @{['other.tar.gz',0]});
	}

	#RAFSeg_chrX
	if(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), @{['rafseg.tar.gz',0]}) == 0){
		_zip_and_tar_fileset($options,
														$RAFSEG_PNG_OUTPUT,
														sprintf($RAFSEG_PNG_TAR,$options->{'tumour_name'}),
														sprintf($RAFSEG_PNG_DIR,$options->{'tumour_name'}),
														$tmp,
														$outdir,
														$contigs,
														$options->{'tumour_name'});
		PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), @{['rafseg.tar.gz',0]});
	}

	#impute input
	if(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), @{['impute_input.tar.gz',0]}) == 0){
		_zip_and_tar_fileset($options,
														$IMPUTE_INPUT_OUTPUT,
														sprintf($IMPUTE_INPUT_TAR,$options->{'tumour_name'}),
														sprintf($IMPUTE_INPUT_DIR,$options->{'tumour_name'}),
														$tmp,
														$outdir,
														undef,
														$options->{'tumour_name'});
		PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), @{['impute_input.tar.gz',0]});
	}

	#impute output
	if(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), @{['impute_output.tar.gz',0]}) == 0){
		_zip_and_tar_fileset($options,
														$IMPUTE_OUTPUT_OUTPUT,
														sprintf($IMPUTE_OUTPUT_TAR,$options->{'tumour_name'}),
														sprintf($IMPUTE_OUTPUT_DIR,$options->{'tumour_name'}),
														$tmp,
														$outdir,
														undef,
														$options->{'tumour_name'});
		PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), @{['impute_output.tar.gz',0]});
	}

	#Cleanup other impute files.
	if(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), @{['cleanup_impute',0]}) == 0){
		for(my $i=1; $i<=$options->{'job_count'}; $i++){
			my $file = File::Spec->catfile($tmp,sprintf($IMPUTE_INPUT,$options->{'tumour_name'},$i));
			next if($options->{'is_male'} && $file =~ m/chr(X|23)/ && ! -e $file); #Skip if male and X results aren't present, this is allowed...
	                foreach my $file_to_del (glob $file){
	                       unlink($file_to_del) or print "Could not unlink $file_to_del";
	                }
		}
		PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), @{['cleanup_impute',0]});
	}

	#LogR and segmented VCF
	if(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), @{['create_seg_vcf',0]}) == 0){
		#Instead of the below we create a VCF and tar gzip and tabix that, copying those across.
		my $baf_seg = File::Spec->catfile($tmp,sprintf($BAF_SEGMENT_TXT,$options->{'tumour_name'}));
		my $log_r_seg = File::Spec->catfile($tmp,sprintf($LOGR_SEGMENT_TXT,$options->{'tumour_name'}));
		my $vcf_seg_out = File::Spec->catfile($tmp,sprintf($SEGMENT_VCF,$options->{'tumour_name'}));
		_generate_segmented_vcf($options,$baf_seg,$log_r_seg,$vcf_seg_out);
		#Now tar gz and tabix
		_bgzip_tabix_vcf($options,$vcf_seg_out);
		#Move tar.gz and tabix to results folder.
		tmp_to_outdir($tmp, $outdir, $tumour_name,
  	              $SEGMENT_VCF_GZ, $SEGMENT_VCF_TABIX);
		PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), @{['create_seg_vcf',0]});
	}

	#Copy number VCF
	if(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), @{['create_cn_vcf',0]}) == 0){
		#Generate cn vcf
		my $cn_file = File::Spec->catfile($tmp,sprintf($CN_VCF,$options->{'tumour_name'}));
		my $subcl_txt = File::Spec->catfile($tmp,sprintf($SUBCLONES_TXT,$options->{'tumour_name'}));
		_run_copy_number_to_vcf($options,$subcl_txt,$cn_file);
		#gz and tabix cn vcf
		_bgzip_tabix_vcf($options,$cn_file);
		#copy gz and tabix to results
		tmp_to_outdir($tmp, $outdir, $tumour_name,
		              $CN_VCF_GZ, $CN_VCF_TABIX);
		PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), @{['create_cn_vcf',0]});
	}

	if(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), @{['single_file_copy',0]}) == 0){
		#Now copy (and tar gz) some per run files.
		my $tumour_dot = $tumour_name;
		$tumour_dot =~ s/[-]/./g;

		#Tumour png
		my $tumour = File::Spec->catfile($tmp,sprintf($TUMOUR_PNG,$tumour_name));
		my $tumour_copy = File::Spec->catfile($outdir,sprintf($TUMOUR_CORRECTED_PNG,$tumour_name));
                if (-e $tumour) {
                       _copy_file($tumour,$tumour_copy);
                }
		#Normal png
		my $normal = File::Spec->catfile($tmp,sprintf($NORMAL_PNG,$tumour_name));
		my $normal_copy = File::Spec->catfile($outdir,sprintf($NORMAL_CORRECTED_PNG,$tumour_name));
                if (-e $normal) {
                        _copy_file($normal,$normal_copy);
                }
		my $subcl_txt = File::Spec->catfile($tmp,sprintf($SUBCLONES_TXT,$options->{'tumour_name'}));
    gzip $subcl_txt => "$subcl_txt.gz" or die "gzip failed: $GzipError\n";
    tmp_to_outdir($tmp, $outdir, $tumour_name,
                  $SUNRISE_PNG, $COPY_NO_PNG, $NON_ROUNDED_PNG, $ALT_COPY_NO_ROUNDED_PNG,
                  $ALT_NON_ROUNDED_PNG, $SECOND_DISTANCE_PNG, $PROFILE_AVERAGE_PNG,
                  $PROFILE_SUBCLONES_PNG, $RHO_PSI_FILE, $NORMAL_CONTAMINATION_FILE,
                  "$SUBCLONES_TXT.gz", $CELLULARITY_PLOIDY_TXT, $STATUS_TXT, $SEED_TXT);

		PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), @{['single_file_copy',0]});
	}
  return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'),0);
}

sub determine_gender {
  my $options = shift;
  my $gender_loci;
  if(defined $options->{'genderloci'}) {
    $gender_loci = $options->{'genderloci'};
  }
  else {
    my $mod_path = dirname(abs_path($0)).'/../share';
    $mod_path = module_dir('Sanger::CGP::Battenberg::Implement') unless(-e File::Spec->catdir($mod_path, 'gender'));

    my $gender_path = File::Spec->catdir($mod_path, 'gender');
    $gender_loci = File::Spec->catfile($gender_path,'GRCh37d5_Y.loci');
  }

  my $command = _which('alleleCounter');
  $command .= sprintf $ALLELE_COUNT_PARA, $options->{'normbam'}, File::Spec->catfile($options->{'tmp'}, 'normal_gender.tsv'), $gender_loci;
  $command .= '-m '.$options->{'mbq'} if exists $options->{'mbq'} && defined $options->{'mbq'};
  system($command);
  my $norm_gender = _parse_gender_results(File::Spec->catfile($options->{'tmp'}, 'normal_gender.tsv'));
  return $norm_gender;
}

sub _parse_gender_results {
  my $file = shift @_;
  my $gender = 'XX';
  open my $fh, '<', $file;
  while(my $line = <$fh>) {
    next if($line =~ m/^#/);
    chomp $line;
    #CHR	POS	Count_A	Count_C	Count_G	Count_T	Good_depth
    my ($chr, $pos, $a, $c, $g, $t, $depth) = split /\t/, $line;
    # all we really care about is the depth
    if($depth > 5) {
      $gender = 'XY';
      last; # presence of ANY male loci in normal is sufficient, we shouldn't be using this to check for 'matchedness'
    }
  }
  close $fh;
  return $gender;
}

sub tmp_to_outdir {
  my ($tmp, $outdir, $tumour_name, @formats) = @_;
  for my $format(@formats) {
    my $from = File::Spec->catfile($tmp,sprintf($format,$tumour_name));
    my $to = File::Spec->catfile($outdir,sprintf($format,$tumour_name));
    if (-e $from) {
      _copy_file($from,$to);
    }
  }
}

sub battenberg_cleanup{
	# uncoverable subroutine
	my $options = shift;
  my $tmp = $options->{'tmp'};
  my $outdir = $options->{'outdir'};
  move ($options->{'logs'},File::Spec->catdir($outdir,'logs')) || die "Error trying to move logs directory '$options->{logs}' -> '".File::Spec->catdir($outdir,'logs')."': $!";

  #delete the entire tmp directory
	remove_tree ($tmp) if(-e $tmp);
	return;
}

sub _run_copy_number_to_vcf{
	my ($options,$subcl_txt,$cn_file) = @_;
	my $tmp = $options->{'tmp'};
	my $command = "$^X ";
  $command .= _which('battenberg_CN_to_VCF.pl');
	$command .= " -o $cn_file";
  $command .= " -r $options->{reference}";
  $command .= " -i $subcl_txt";
  $command .= " -ra $options->{assembly}";
  $command .= " -rs $options->{species}";
  $command .= " -msq $options->{protocol} -wsq $options->{protocol}";
  $command .= " -msp $options->{platform} -wsp $options->{platform}";
  if(defined $options->{'allele-counts'}) {
    # these are sample names when allele-count is in effect
    $command .= " -msn $options->{tumbam}";
    $command .= " -wsn $options->{normbam}";
  }
  else {
    $command .= " -sbm $options->{tumbam}";
    $command .= " -sbw $options->{normbam}";
  }
	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
  return;
}

sub _bgzip_tabix_vcf{
	my ($options,$vcf) = @_;
	my $tmp = $options->{'tmp'};
	my $vcf_gz = $vcf.'.gz';
  my $bgzip = _which('bgzip');
  $bgzip .= sprintf ' -c %s > %s', $vcf, $vcf_gz;

  my $tabix = _which('tabix');
  $tabix .= sprintf ' -p vcf %s', $vcf_gz;

  my @commands = [$bgzip, $tabix];

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $bgzip, 0);
  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $tabix, 0);
  return;
}

sub _generate_segmented_vcf{
	my ($options,$baf,$logr,$vcf_out) = @_;
	my $stripped_baf = basename($baf);
	my $stripped_logr = basename($logr);

	open (my $fh_BAF, '<',$baf);
	open (my $fh_logR,'<',$logr);
	open (my $fhw, '>', $vcf_out);
		my $vcf=Vcf->new();
		$vcf->add_header_line({key=>'INFO',ID=>'BAFP', Number=>'1',Type=>'Float',Description=>"BAF:BAF values after phasing with Impute"});
		$vcf->add_header_line({key=>'INFO',ID=>'BAFP', Number=>'1',Type=>'Float',Description=>"BAF:BAF values after phasing with Impute"});
		$vcf->add_header_line({key=>'INFO',ID=>'BAFE', Number=>'1',Type=>'Float',Description=>"BAFphased:BAF values after correction for errors in phasing"});
		$vcf->add_header_line({key=>'INFO',ID=>'BAFS', Number=>'1',Type=>'Float',Description=>"BAFseg:segmented BAF values for each heterozygous SNP"});
		$vcf->add_header_line({key=>'INFO',ID=>'LOGR', Number=>'1',Type=>'Float',Description=>"logRseg:segmented logR values for each heterozygous SNP"});
		$vcf->add_header_line({key=>'SAMPLE',ID=>'BAMNAME', Number=>'0',Type=>'String',SampleName=>"$options->{tumour_name}"});
		$vcf->add_header_line({key=>'Data',value=>"combined data from $stripped_baf and $stripped_logr"});
		$vcf->add_columns($options->{'tumour_name'});

		print $fhw $vcf->format_header();
		my %logR;
		while(<$fh_logR>) {
			my @f=split(/\s+/,$_);
      #use chr + pos as a unique index
      my $index = join ':', $f[0], $f[1];
		 	$logR{$index}=$f[2];
		}
		close($fh_logR);
		while (<$fh_BAF>) {
			chomp;
			my @f=split(' ',$_);
      next if ($f[0] =~ /^Chromosome/);
      my %out;
      #use chr + pos as a unique index
      my $index = join ':', $f[0], $f[1];
     if (!exists $logR{$index}){
        $logR{$index}="NA";
      }
      $out{CHROM}  = $f[0];
      $out{POS}    = $f[1];
      $out{ID}     = '.';
      $out{ALT}    = ['.'];
      $out{REF}    = '.';
      $out{QUAL}   = '.';
      $out{FILTER} = ['.'];
      $out{INFO}   = { BAFP=>$f[2], BAFE=>$f[3], BAFS=>$f[4], LOGR=>$logR{$index} };
      print $fhw $vcf->format_line(\%out) or die("Error printing line to vcf file.");
		}

		$vcf->close();

	close($fh_BAF);
	close($fhw);
	return;
}

sub _copy_file{
	my ($from,$to) = @_;
	PCAP::Cli::file_for_reading('From location copy check',$from);
	die("Error: Failed to copy file $from to $to.") if(!copy($from,$to));
	return;
}

sub _zip_and_tar_fileset{
	my ($options,$filepattern,$tar,$dir,$location,$copyloc,$file_match_list,$sample_name) = @_;
	my @files = ();
	#Get a list of files matching the pattern
	for(my $i=0; $i<$options->{'job_count'}; $i++){
		my $file;
		if(defined($file_match_list)){
			$file = File::Spec->catfile($location,sprintf($filepattern,$sample_name,$file_match_list->[$i]));
			next if($options->{'is_male'} && $file =~ m/chr(X|23)/ && ! -e $file); #Skip if male and X results aren't present, this is allowed...
		}else{
			$file = File::Spec->catfile($location,sprintf($filepattern,$sample_name,$i+1));
			next if($options->{'is_male'} && $file =~ m/chr(X|23)/ && ! -e $file); #Skip if male and X results aren't present, this is allowed...
		}
		if (-e $file) {
		  PCAP::Cli::file_for_reading('file for tar.gz',$file);
		  push(@files,$file);
		}
	}

	my $tarball = _targzFileSet($options,\@files,$tar,$dir);
	my $name = fileparse($tarball);
	my $copied_loc = File::Spec->catfile($copyloc,$name);
	die("Error: Failed to copy file $tarball to $copied_loc.") if(!copy($tarball,$copied_loc));
	push(@files,$tarball);
	foreach my $file_to_del (@files){
		unlink($file_to_del) or print "Could not unlink $file_to_del";
	}
	return;
}

sub _targzFileSet{
	my ($options,$fileList,$archive_name,$folder_name) = @_;
	my $tmp = $options->{'tmp'};
	my $dir = File::Spec->catdir($options->{'tmp'},$folder_name);
	if(! -e $dir){
		mkdir($dir) or die("Error trying to make directory $dir\n");
	}
	#Iterate through each file and copy into the folder.
	foreach my $file_to_cp(@$fileList){
	  my $fname = fileparse($file_to_cp);
		my $copied_loc = File::Spec->catfile($dir,$fname);
		die("Error: Failed to copy file $file_to_cp to $copied_loc.") if(!copy($file_to_cp,$copied_loc));
	}

	my $command = "cd $tmp; ".q{tar -czf }.$archive_name;
	$command .= q{ -C }.$tmp;
	$command .= q{ }.$folder_name;

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
	return File::Spec->catfile($options->{'tmp'},$archive_name);
}

sub _writeNormalContaminationToFile{
	# uncoverable subroutine
	my $options = shift;
	my $outdir = $options->{'outdir'};
  my $tmp = $options->{'tmp'};
	my $sample_name = $options->{'tumour_name'};
	my $normal_cont = $options->{'normc'};
	my $file = File::Spec->catfile($tmp,sprintf($NORMAL_CONTAMINATION_FILE,$sample_name));
	my $FH;
	open($FH, '>', $file) || die("Error opening file '$file' for write: $!");
		print $FH ($normal_cont,"\n");
	close($FH);
	return;
}

sub _writeSeedToFile{
  # uncoverable subroutine
  my $options = shift;
  my $outdir = $options->{'outdir'};
  my $tmp = $options->{'tmp'};
  my $sample_name = $options->{'tumour_name'};
  my $seed = $options->{'seed'};
  my $file = File::Spec->catfile($tmp,sprintf($SEED_TXT,$sample_name));
  my $FH;
  open($FH, '>', $file) || die("Error opening file '$file' for write: $!");
  print $FH ($seed,"\n");
  close($FH);
  return;
}

sub _extractRhoPsiFromFile{
	# uncoverable subroutine
	my $options = shift;
	my $outdir = $options->{'outdir'};
  my $tmp = $options->{'tmp'};
	my $sample_name = $options->{'tumour_name'};
  my ($name,$rho,$psi);
  my $file_path = File::Spec->catfile($tmp,sprintf($RHO_PSI_FILE,$sample_name));
	my $FH;
	open($FH,'<',$file_path) or die("Error trying to read rho_psi_file '$file_path': $!");
		while(<$FH>){
			my $line = $_;
			next unless($line =~ m/^FRAC_GENOME/);
			chomp($line);
			($name,$rho,$psi,undef) = split(/\t/,$line);
		}
	close($FH) or die("Error trying to close rho_psi_file '$file_path': $!");
	return ($rho,$psi);

}

sub file_line_count_with_ignore{
	my ($file,$ignore_contigs) = @_;
	my $contig_count = 0;
  open my $FH, '<', $file or die("Error trying to open $file: $!\n");
  while(<$FH>){
    my $line = $_;
    chomp($line);
    my ($contig,undef) = split(/\s+/,$line);
    next if(first {$_ eq $contig} @{$ignore_contigs});
    $contig_count++;
  }
  close($FH);
  return $contig_count;
}

sub read_contigs_from_file_with_ignore{
	my ($file,$ignore_contigs) = @_;
	my @contigs = ();
	{
		my $FH;
		open($FH, '<', $file) or die("Error trying to open file: $!\n");
			while(<$FH>){
    		my $line = $_;
    		chomp($line);
    		my ($con,undef) = split(/\s+/,$line);
    		$con =~ s/chr//;  # handle hg19, removing chr prefix
    		my $match=0;
    		foreach my $ign(@$ignore_contigs){
    			if("$ign" eq "$con"){
    				$match = 1;
    				last;
    			}
    		}
    		next if($match);
    		push(@contigs,$con);
    	}
		close($FH);
	}
	return \@contigs;
}

sub read_contigs_from_file{
	my $file = shift;
	my @contigs = ();
	{
		my $FH;
		open($FH, '<', $file) or die("Error trying to open file: $!\n");
			while(<$FH>){
    		my $line = $_;
    		chomp($line);
    		my ($con,undef) = split(/\s+/,$line);
    		push(@contigs,$con);
    	}
		close($FH);
	}
	return \@contigs;
}

sub _which {
  my $prog = shift;
  my $l_bin = $Bin;
  my $path = File::Spec->catfile($l_bin, $prog);
  $path = which($prog) unless(-e $path);
  die "Failed to find $prog in PATH or local bin folder" unless(defined $path);
  return $path;
}

sub limited_indices {
  my ($options, $index_in, $count) = @_;
  my @indices;
  if(exists $options->{'limit'}) {
    # main script checks index is not greater than limit or < 1
    my $base = $index_in;
    while($base <= $count) {
      push @indices, $base;
      $base += $options->{'limit'};
    }
  }
  else {
    push @indices, $index_in;
  }
  return @indices;
}

sub _find_num_loci_per_file {
  my ($thousand_genome_loci_dir, $ignored_contigs) = @_;

  #Find the number of loci (lines) in each 1000 genome file
  my $loci_pattern = File::Spec->catfile($thousand_genome_loci_dir, $ONEKGEN_LOCI_FILE_PATTERN);

  #Convert ignored_contigs array into a hash
#  my %ignored_contigs_hash = map { $_ => 1 } @$ignored_contigs;
  my %ignored_contigs_hash;
  foreach my $contig (@$ignored_contigs) {
    $contig =~ s/^chr//;
    $ignored_contigs_hash{$contig} = 1;
  }

  my @loci_files = < $loci_pattern >;

  my %lines_per_file;
  my $total = 0;
  foreach my $file (@loci_files) {
    open my $fh, $file or die "Unable to open $file";
    #check the first line for the actual chromosome name
    my $first_line = <$fh>;
    my ($chr, $pos) = split " ", $first_line;
    $chr =~ s/^chr//;

    #skip if we want to ignore this file
    next if ($ignored_contigs_hash{$chr});

    while ( <$fh> ) {};
    my $count = $.;
    $lines_per_file{$file} = $count;
    $total += $count;
    close $fh;
  }

  #Return the number of loci (lines) for each file and the grand total over all files
  return (\%lines_per_file, $total);
}

sub _create_split_files {
  my ($file, $total_loci, $number_of_split_files, $tmpdir, $first_number) = @_;

  #Maximum number of lines to write in each file
  my $number_of_split_loci = ceil($total_loci / $number_of_split_files);

  my $line_count = 0;
  my $file_count = $first_number + 1; #start counting at 1 rather than 0

  #Name the split files after the parent file
  my ($basename, $path, $suffix) = fileparse($file, qr/\.[^.]*/);
  open my $FH, $file or die "Unable to open $file";

  my $split_filename = $basename . "_split$file_count" . $suffix;
  my $split_file = File::Spec->catfile($tmpdir, $split_filename);
  open my $SPLIT_FH, '>', $split_file or die "Unable to open $split_file";

  while (my $line = <$FH>) {
    #Write the correct number of lines to the split file
    if ($line_count >= $number_of_split_loci) {
      close $SPLIT_FH;
      $file_count++;
      $split_file = File::Spec->catfile($tmpdir, $basename . "_split$file_count" . $suffix);
      open $SPLIT_FH, '>', $split_file or die "Unable to open $split_file";
      $line_count = 0;
    }
    print $SPLIT_FH $line;
    $line_count++;
  }

  close $SPLIT_FH;
  close $FH;
}

sub _lociNameMap {
  my ($kgenloc) = @_;

  my $loci_pattern = File::Spec->catfile($kgenloc, $ONEKGEN_LOCI_FILE_PATTERN);
  my @loci_files = < $loci_pattern >;

  my $loci_names_to_index;
  foreach my $file (@loci_files) {
    open my $fh, $file or die "Unable to open $file";
    my ($loci_index) = $file =~ /$ONEKGEN_LOCI_FILE_REGEX/;
    #check the first line for the actual chromosome name
    my $first_line = <$fh>;
    my ($chr, $pos) = split " ", $first_line;
    $chr =~ s/^chr//;
    $loci_names_to_index->{$chr} = $loci_index;

    close $fh;
  }
  return $loci_names_to_index;
}

sub _get_chrom_names {
  my ($impute_info_file, $is_male) = @_;

  my @chroms;
  open my $impute_fh, $impute_info_file or die "Unable to open $impute_info_file";
  while (my $line = <$impute_fh>) {
    chomp $line;
    my ($chrom, $impute_legend, $genetic_map, $impute_hap, $start, $end, $is_par) = split "\t", $line;
    if ($is_male && $is_par) {
      push @chroms, $chrom;
    } else {
      push @chroms, $chrom;
    }
  }
  close $impute_fh;
  my @unique_chroms = uniq(@chroms);
  return \@unique_chroms;
}

sub _get_chroms_as_string {
  my ($impute_info_file, $is_male) = @_;

  my $chroms = _get_chrom_names($impute_info_file, $is_male);
  my @chroms_string = map {'"' . $_ . '"'} @$chroms;
  my $chr_names = join ',', @chroms_string;
  return $chr_names;
}

1;
