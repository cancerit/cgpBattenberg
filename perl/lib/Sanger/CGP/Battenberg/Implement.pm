package Sanger::CGP::Battenberg::Implement;

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
use Vcf;

use File::ShareDir qw(module_dir);

use Sanger::CGP::Battenberg;

use PCAP::Threaded;
use PCAP::Bam;

use Data::Dumper;

const my $ALLELE_COUNT_CMD => q{ -l %s -b %s -o %s -m %d };
const my $RUN_BAF_LOG => q{ %s/RunBAFLogR.R %s %s %s_alleleFrequencies.txt %s_alleleFrequencies.txt %s_ 10 %s %s };
const my $IMPUTE_FROM_AF => q{ %s/GenerateImputeInputFromAlleleFrequencies.R %s %s %s %s_alleleFrequencies.txt %s_alleleFrequencies.txt %s_impute_input_chr %d %s};
const my $RUN_IMPUTE => q{ %s/RunImpute.R %s %s %s %s %s_impute_input_chr %s_impute_output_chr %d};
const my $COMBINE_IMPUTE => q{ %s/CombineImputeOutputs.R %s %s %s %s_impute_output_chr %d};
const my $HAPLOTYPE_BAF => q{ %s/RunGetHaplotypedBAFs.R %s %s %s %s_alleleFrequencies.txt %s_alleleFrequencies.txt %s_impute_output_chr %s %s_ %d};
const my $PLOT_HAPLOTYPE_BAFS => q{ %s/PlotHaplotypedData.R %s %s %d %s_chr%d_heterozygousMutBAFs_haplotyped.txt %s};
const my $COMBINE_BAFS => q{ %s/CombineBAFfiles.R %s %s %s %s_chr _heterozygousMutBAFs_haplotyped.txt %s_allChromosomes_heterozygousMutBAFs_haplotyped.txt};
const my $SEGMENT_PHASED => q{ %s/segmentBAFphased.R %s %s %d %d};
const my $FIT_COPY_NUMBER => q{ %s/FitCopyNumber.R %s %s %s_ %d %d %f %d %f %f %f %f};
const my $CALL_SUBCLONES => q{ %s/callSubclones.R %s %s %s %s %d %d};

const my $ALLELE_COUNT_OUTPUT => q{%s_alleleFrequencies_chr%d.txt};
const my $ALLELE_LOCI_NAME => q{1000genomesloci2012_chr%d.txt};
const my $ALLELE_COUNT_TAR => q{%s_allelecounts.tar.gz};
const my $ALLELE_COUNT_DIR => q{%s_allelecounts};
const my $ALLELE_COUNT_SCRIPT => q{alleleCounter};
const my $NORMAL_CONTAMINATION_FILE => q{%s_normal_contamination.txt};
const my $RHO_PSI_FILE => q{%s_rho_and_psi.txt};
const my $SUBCLONE_PNG_OUTPUT => q{%s_subclones_chr%s.png};
const my $SUBCLONE_TAR => q{%s_subclones.tar.gz};
const my $SUBCLONE_DIR => q{%s_subclones};
const my $HETDATA_OUTPUT => q{%s_chr%s_heterozygousData.png};
const my $HETDATA_TAR => q{%s_hetdata.tar.gz};
const my $HETDATA_DIR => q{%s_hetdata};
const my $HETBAFTXT_OUTPUT => q{%s_chr%s_heterozygousMutBAFs_haplotyped.txt};
const my $HETBAFTXT_TAR => q{%s_hetbaf.tar.gz};
const my $HETBAFTXT_DIR => q{%s_hetbaf};
const my $OTHER_PNG_OUTPUT => q{%s_chr%s.png};
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
const my $TUMOUR_PNG	=> q{%s_Tumor%s.png};
const my $NORMAL_PNG	=> q{%s_Germline%s.png};
const my $COPY_NO_PNG	=> q{%s_second_copynumberprofile.png};
const my $NON_ROUNDED_PNG	=> q{%s_second_nonroundedprofile.png};
const my $ALT_COPY_NO_ROUNDED_PNG => q{%s_copynumberprofile.png};
const my $ALT_NON_ROUNDED_PNG => q{%s_nonroundedprofile.png};
const my $SECOND_DISTANCE_PNG => q{%s_second_distance.png};
const my $BAF_SEGMENT_TXT => q{%s.BAFsegmented.txt};
const my $LOGR_SEGMENT_TXT => q{%s.logRsegmented.txt};
const my $SUBCLONES_TXT => q{%s_subclones.txt};
const my $SEGMENT_VCF => q{%s_logR_Baf_segmented.vcf};
const my $SEGMENT_VCF_GZ => q{%s_logR_Baf_segmented.vcf.gz};
const my $SEGMENT_VCF_TABIX => q{%s_logR_Baf_segmented.vcf.gz.tbi};
const my $CN_VCF => q{%s_battenberg_cn.vcf};
const my $CN_VCF_GZ => q{%s_battenberg_cn.vcf.gz};
const my $CN_VCF_TABIX => q{%s_battenberg_cn.vcf.gz.tbi};

const my $RSCRIPT => q{Rscript};
const my $IMPUTE_EXE => q{impute2};

const my $IMPUTE_INPUT => q{%s_impute_input_chr%d_*K.*};
const my $IMPUTE_OUTPUT => q{%s_impute_output_chr%d.txt};

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

sub prepare {
  my $options = shift;
  $options->{'tumbam'} = File::Spec->rel2abs($options->{'tumbam'});
	$options->{'normbam'} = File::Spec->rel2abs($options->{'normbam'});
  $options->{'tumour_name'} = (PCAP::Bam::sample_name($options->{'tumbam'}))[0];
  $options->{'normal_name'} = (PCAP::Bam::sample_name($options->{'normbam'}))[0];
  my $mod_path = dirname(abs_path($0)).'/../share';
  $mod_path = module_dir('Sanger::CGP::Battenberg::Implement') unless(-e File::Spec->catdir($mod_path, 'battenberg'));
	$options->{'mod_path'} = $mod_path;
	$options->{'bat_path'} = File::Spec->catdir($mod_path, 'battenberg');
	$options->{'tmp'} = File::Spec->rel2abs($options->{'tmp'});
  return 1;
}

sub battenberg_allelecount{
	# uncoverable subroutine
	my ($index, $options) = @_;
	return 1 if(exists $options->{'index'} && $index != $options->{'index'});
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);

	my $k_gen_loc = File::Spec->rel2abs($options->{'1kgenloc'});

	#A little sorting so we get all allele counting done in parallel (we have 2 * no of contigs threads)
  my $input = $options->{'tumbam'};
  my $sname = $options->{'tumour_name'};
	my $lookup = $index;

	if($index>$options->{'job_count'}){
		$lookup = $index - $options->{'job_count'};
		$input = $options->{'normbam'};
		$sname = $options->{'normal_name'};
	}

	my $loci_file = File::Spec->catfile($k_gen_loc,sprintf($ALLELE_LOCI_NAME,$lookup));
	my $alleleCountOut = File::Spec->rel2abs(File::Spec->catfile($tmp,sprintf($ALLELE_COUNT_OUTPUT,$sname,$lookup)));

	my $command = "cd $tmp; "._which($ALLELE_COUNT_SCRIPT) || die "Unable to find $ALLELE_COUNT_SCRIPT in path";

	$command .= sprintf($ALLELE_COUNT_CMD,
							$loci_file,
							$input,
							$alleleCountOut,
							$options->{'mbq'}
							);

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);

	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);

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

  $command .= sprintf($RUN_BAF_LOG,
  								$options->{'bat_path'}, #Used to build path to run Rscript
  								$options->{'bat_path'}, #Used by rScript to path to modules
									$impute_info,
									$options->{'tumour_name'},
									$options->{'normal_name'},
									$options->{'tumour_name'},
									$options->{'tumour_name'},
									$thou_gen_loc,
  						);
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
	#Use the share directory so we don't have to symlink anything.
	my $command = "cd $tmp; ";
	$command .= _which($RSCRIPT) || die "Unable to find $RSCRIPT in path";

 	$command .= sprintf($IMPUTE_FROM_AF,
												$options->{'bat_path'},
												$options->{'bat_path'},
												$impute_info,
												$options->{'is_male'},
												$options->{'tumour_name'},
												$options->{'normal_name'},
												$options->{'tumour_name'},
												$index,
												$prob_loci
 											);

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

	#Use the share directory so we don't have to symlink anything.
	my $command = "cd $tmp; ";
	$command .= _which($RSCRIPT) || die "Unable to find $RSCRIPT in path";

	$command .= sprintf($RUN_IMPUTE,
												$options->{'bat_path'},
												$options->{'bat_path'},
												$impute_info,
												$impute_exe,
												$options->{'is_male'},
												$options->{'tumour_name'},
												$options->{'tumour_name'},
												$index,
											);

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

	#Use the share directory so we don't have to symlink anything.
	my $command = "cd $tmp; ";
	$command .= _which($RSCRIPT) || die "Unable to find $RSCRIPT in path";

	$command .= sprintf($COMBINE_IMPUTE,
												$options->{'bat_path'},
												$options->{'bat_path'},
												$impute_info,
												$options->{'is_male'},
												$options->{'tumour_name'},
												$index,
												);

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

	my $command = "cd $tmp; ";
	$command .= _which($RSCRIPT) || die "Unable to find $RSCRIPT in path";

	$command .= sprintf($HAPLOTYPE_BAF,
												$options->{'bat_path'},
												$options->{'bat_path'},
												$impute_info,
												$options->{'is_male'},
												$options->{'tumour_name'},
												$options->{'normal_name'},
												$options->{'tumour_name'},
												$options->{'tumour_name'},
												$options->{'tumour_name'},
												$index,
							);

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

	my $command = "cd $tmp; ";
	$command .= _which($RSCRIPT) || die "Unable to find $RSCRIPT in path";

	$command .= sprintf($PLOT_HAPLOTYPE_BAFS,
												$options->{'bat_path'},
												$options->{'bat_path'},
												$impute_info,
												$index,
												$options->{'tumour_name'},
												$index,
												$options->{'tumour_name'},
											);

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
	my $command = "cd $tmp; ";
	$command .= _which($RSCRIPT) || die "Unable to find $RSCRIPT in path";

 	$command .= sprintf($COMBINE_BAFS,
												$options->{'bat_path'},
												$options->{'bat_path'},
												$impute_info,
												$options->{'is_male'},
												$options->{'tumour_name'},
												$options->{'tumour_name'},
 								);

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub battenberg_segmentphased{
	# uncoverable subroutine
	my $options = shift;
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
	my $command = "cd $tmp; ";
	$command .= _which($RSCRIPT) || die "Unable to find $RSCRIPT in path";

 	$command .= sprintf($SEGMENT_PHASED,
												$options->{'bat_path'},
												$options->{'bat_path'},
												$options->{'tumour_name'},
												$options->{'seg_gamma'},
												$options->{'phase_gamma'},
										);

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);

}

sub battenberg_fitcopyno{
	# uncoverable subroutine
	my $options = shift;
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
	my $command = "cd $tmp; ";
	$command .= _which($RSCRIPT) || die "Unable to find $RSCRIPT in path";
	$command .= sprintf($FIT_COPY_NUMBER,
												$options->{'bat_path'},
												$options->{'bat_path'},
												$options->{'tumour_name'},
												$options->{'tumour_name'},
												$options->{'clonality_dist'},
												$options->{'ascat_dist'},
												$options->{'balanced_thresh'},
												$options->{'plat_gamma'},
												$options->{'min_ploidy'},
												$options->{'max_ploidy'},
												$options->{'min_rho'},
												$options->{'min_goodness'},
											);

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub battenberg_callsubclones{
	# uncoverable subroutine
	my $options = shift;
	my $tmp = $options->{'tmp'};
	my $impute_info = $options->{'impute_info'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
	my $command = "cd $tmp; ";
	$command .= _which($RSCRIPT) || die "Unable to find $RSCRIPT in path";

	$command .= sprintf($CALL_SUBCLONES,
												$options->{'bat_path'},
												$options->{'bat_path'},
												$impute_info,
												$options->{'is_male'},
												$options->{'tumour_name'},
												$options->{'seg_gamma'},
												$options->{'plat_gamma'},
									);

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub battenberg_finalise{
	# uncoverable subroutine
	my $options = shift;
	my $outdir = $options->{'outdir'};
  my $tmp = $options->{'tmp'};

  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

	my ($rho,$psi) = _extractRhoPsiFromFile($options);

	#Calculate normal contamination
	my $normal_contamination = 2*(1-$rho) /
												(2*(1-$rho) + $psi * $rho);
	$options->{'normc'} = $normal_contamination;
	_writeNormalContaminationToFile($options);

	#Tarball allelcounts and copy to results folder
	if (PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), @{['allele_tar_gz',0]}) == 0){
		_zip_and_tar_fileset($options,
														$ALLELE_COUNT_OUTPUT,
														sprintf($ALLELE_COUNT_TAR,$options->{'normal_name'}),
														sprintf($ALLELE_COUNT_DIR,$options->{'normal_name'}),
														$tmp,
														$outdir,
														undef);

		_zip_and_tar_fileset($options,
														$ALLELE_COUNT_OUTPUT,
														sprintf($ALLELE_COUNT_TAR,$options->{'tumour_name'}),
														sprintf($ALLELE_COUNT_DIR,$options->{'tumour_name'}),
														$tmp,
														$outdir,
														undef);

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
														$contigs);
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
														$contigs);
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
														undef);
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
														$contigs);
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
														$contigs);
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
														undef);
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
														undef);
		PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), @{['impute_output.tar.gz',0]});
	}

	#Cleanup other impute files.
	if(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), @{['cleanup_impute',0]}) == 0){
		for(my $i=1; $i<=$options->{'job_count'}; $i++){
			my $file = File::Spec->catfile($tmp,sprintf($IMPUTE_INPUT,$options->{'tumour_name'},$i));
			unlink glob $file;
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
		my $gz_file = File::Spec->catfile($tmp,sprintf($SEGMENT_VCF_GZ,$options->{'tumour_name'}));
		my $gz_file_copy = File::Spec->catfile($outdir,sprintf($SEGMENT_VCF_GZ,$options->{'tumour_name'}));
		my $tabix_file = File::Spec->catfile($tmp,sprintf($SEGMENT_VCF_TABIX,$options->{'tumour_name'}));
		my $tabix_file_copy = File::Spec->catfile($outdir,sprintf($SEGMENT_VCF_TABIX,$options->{'tumour_name'}));
		_copy_file($gz_file,$gz_file_copy);
		_copy_file($tabix_file,$tabix_file_copy);
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
		my $gz_file = File::Spec->catfile($tmp,sprintf($CN_VCF_GZ,$options->{'tumour_name'}));
		my $gz_file_copy = File::Spec->catfile($outdir,sprintf($CN_VCF_GZ,$options->{'tumour_name'}));
		my $tabix_file = File::Spec->catfile($tmp,sprintf($CN_VCF_TABIX,$options->{'tumour_name'}));
		my $tabix_file_copy = File::Spec->catfile($outdir,sprintf($CN_VCF_TABIX,$options->{'tumour_name'}));
		_copy_file($gz_file,$gz_file_copy);
		_copy_file($tabix_file,$tabix_file_copy);
		PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), @{['create_cn_vcf',0]});
	}

	if(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), @{['single_file_copy',0]}) == 0){
		#Now copy (and tar gz) some per run files.
		#Sunrise plot
		my $sunrise = File::Spec->catfile($tmp,sprintf($SUNRISE_PNG,$options->{'tumour_name'}));
		my $sunrise_copy = File::Spec->catfile($outdir,sprintf($SUNRISE_PNG,$options->{'tumour_name'}));
		_copy_file($sunrise,$sunrise_copy);
		#Tumour png
		my $tumour = File::Spec->catfile($tmp,sprintf($TUMOUR_PNG,$options->{'tumour_name'},$options->{'tumour_name'}));
		my $tumour_copy = File::Spec->catfile($outdir,sprintf($TUMOUR_PNG,$options->{'tumour_name'},$options->{'tumour_name'}));
		_copy_file($tumour,$tumour_copy);
		#Normal png
		my $normal = File::Spec->catfile($tmp,sprintf($NORMAL_PNG,$options->{'tumour_name'},$options->{'tumour_name'}));
		my $normal_copy = File::Spec->catfile($outdir,sprintf($NORMAL_PNG,$options->{'tumour_name'},$options->{'tumour_name'}));
		_copy_file($normal,$normal_copy);
		#Copy no
		my $copyno = File::Spec->catfile($tmp,sprintf($COPY_NO_PNG,$options->{'tumour_name'}));
		my $copyno_copy = File::Spec->catfile($outdir,sprintf($COPY_NO_PNG,$options->{'tumour_name'}));
		_copy_file($copyno,$copyno_copy);
		#Non rounded copy number
		my $nonrounded = File::Spec->catfile($tmp,sprintf($NON_ROUNDED_PNG,$options->{'tumour_name'}));
		my $nonrounded_copy = File::Spec->catfile($outdir,sprintf($NON_ROUNDED_PNG,$options->{'tumour_name'}));
		_copy_file($nonrounded,$nonrounded_copy);
		#alt_non_rounded copy number
		my $alt_non_rounded = File::Spec->catfile($tmp,sprintf($ALT_COPY_NO_ROUNDED_PNG,$options->{'tumour_name'}));
		my $alt_non_rounded_copy = File::Spec->catfile($outdir,sprintf($ALT_COPY_NO_ROUNDED_PNG,$options->{'tumour_name'}));
		_copy_file($alt_non_rounded,$alt_non_rounded_copy);
		#Alt cn
		my $alt_cn = File::Spec->catfile($tmp,sprintf($ALT_NON_ROUNDED_PNG,$options->{'tumour_name'}));
		my $alt_cn_copy = File::Spec->catfile($outdir,sprintf($ALT_NON_ROUNDED_PNG,$options->{'tumour_name'}));
		_copy_file($alt_cn,$alt_cn_copy);
		#Second distance
		my $sec_dist = File::Spec->catfile($tmp,sprintf($SECOND_DISTANCE_PNG,$options->{'tumour_name'}));
		my $sec_dist_copy = File::Spec->catfile($outdir,sprintf($SECOND_DISTANCE_PNG,$options->{'tumour_name'}));
		_copy_file($sec_dist,$sec_dist_copy);
		#Rho psi txt
		my $rho_psi = File::Spec->catfile($tmp,sprintf($RHO_PSI_FILE,$options->{'tumour_name'}));
		my $rho_psi_copy = File::Spec->catfile($outdir,sprintf($RHO_PSI_FILE,$options->{'tumour_name'}));
		_copy_file($rho_psi,$rho_psi_copy);
		#normal contamination txt
		my $norm_cont = File::Spec->catfile($tmp,sprintf($NORMAL_CONTAMINATION_FILE,$options->{'tumour_name'}));
		my $norm_cont_copy = File::Spec->catfile($outdir,sprintf($NORMAL_CONTAMINATION_FILE,$options->{'tumour_name'}));
		_copy_file($norm_cont,$norm_cont_copy);
		#subclones txt
		my $subcl_txt = File::Spec->catfile($tmp,sprintf($SUBCLONES_TXT,$options->{'tumour_name'}));
		my $subcl_txt_copy = File::Spec->catfile($outdir,sprintf($SUBCLONES_TXT,$options->{'tumour_name'}));
		_copy_file($subcl_txt,$subcl_txt_copy);
		PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), @{['single_file_copy',0]});
	}
	#Lastly move the logs file
	if(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), @{['move_log_dir',0]}) == 0){
			move ($options->{'logs'},File::Spec->catdir($outdir,'logs'))
      	|| die "Error trying to move logs directory '$options->{logs}' -> '".File::Spec->catdir($outdir,'logs')."': $!";
		PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), @{['move_log_dir',0]});
	}
  return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'),0);
}

sub battenberg_cleanup{
	# uncoverable subroutine
	my $options = shift;
  my $tmp = $options->{'tmp'};
  #delete the entire tmp directory
	remove_tree ($tmp);
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
  $command .= " -sbm $options->{tumbam}";
  $command .= " -sbw $options->{normbam}";
  $command .= " -ra $options->{assembly}" if(defined $options->{'assembly'});
  $command .= " -rs $options->{species}" if(defined $options->{'species'});
  $command .= " -msq $options->{protocol} -wsq $options->{protocol}" if(defined $options->{'protocol'});
  $command .= " -msp $options->{platform} -wsp $options->{platform}" if(defined $options->{'platform'});

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

  my @commands = ($bgzip, $tabix);

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), \@commands, 0);
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
		 	$logR{$f[0]}=$f[3];
		}
		close($fh_logR);
		while (<$fh_BAF>) {
			chomp;
			my @f=split(' ',$_);
			if ($f[0]=~/snp/){
				my %out;
				if (!exists $logR{$f[0]}){
					$logR{$f[0]}="NA";
				}
				$out{CHROM}  = $f[1];
				$out{POS}    = $f[2];
				$out{ID}     = $f[0];
				$out{ALT}    = ['.'];
				$out{REF}    = '.';
				$out{QUAL}   = '.';
				$out{FILTER} = ['.'];
				$out{INFO}   = { BAFP=>$f[3], BAFE=>$f[4], BAFS=>$f[5], LOGR=>$logR{$f[0]} };
				print $fhw $vcf->format_line(\%out) or die("Error printing line to vcf file.");
			}
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
	my ($options,$filepattern,$tar,$dir,$location,$copyloc,$file_match_list) = @_;
	my @files = ();
	#Get a list of files matching the pattern
	for(my $i=0; $i<$options->{'job_count'}; $i++){
		my $file;
		if(defined($file_match_list)){
			$file = File::Spec->catfile($location,sprintf($filepattern,$options->{'tumour_name'},$file_match_list->[$i]));
		}else{
			$file = File::Spec->catfile($location,sprintf($filepattern,$options->{'tumour_name'},$i+1));
		}
		PCAP::Cli::file_for_reading('file for tar.gz',$file);
		push(@files,$file);
	}

	my $tarball = _targzFileSet($options,\@files,$tar,$dir);
	my $name = fileparse($tarball);
	my $copied_loc = File::Spec->catfile($copyloc,$name);
	die("Error: Failed to copy file $tarball to $copied_loc.") if(!copy($tarball,$copied_loc));
	push(@files,$tarball);
	foreach my $file_to_del (@files){
		unlink($file_to_del);
	}
	return;
}

sub _targzFileSet{
	my ($options,$fileList,$archive_name,$folder_name) = @_;
	my $tmp = $options->{'tmp'};
	my $dir = File::Spec->catdir($options->{'tmp'},$folder_name);
	if(! -e $dir){
		mkdir($dir) or croak("Error trying to make directory $dir\n");
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

sub _extractRhoPsiFromFile{
	# uncoverable subroutine
	my $options = shift;
	my $outdir = $options->{'outdir'};
  my $tmp = $options->{'tmp'};
	my $sample_name = $options->{'tumour_name'};
  my ($name,$rho,$psi);
  my $file_path = File::Spec->catfile($tmp,sprintf($RHO_PSI_FILE,$sample_name));
	my $FH;
	open($FH,'<',$file_path) or croak("Error trying to read rho_psi_file '$file_path': $!");
		while(<$FH>){
			my $line = $_;
			next unless($line =~ m/^FRAC_GENOME/);
			chomp($line);
			($name,$rho,$psi,undef) = split(/\t/,$line);
		}
	close($FH) or croak("Error trying to close rho_psi_file '$file_path': $!");
	return ($rho,$psi);

}

sub file_line_count_with_ignore{
	my ($file,$ignore_contigs) = @_;
	my $contig_count = 0;
  {
    my $FH;
    open($FH, '<', $file) or die("Error trying to open $file: $!\n");
    	while(<$FH>){
    		my $line = $_;
    		chomp($line);
    		my ($contig,undef) = split(/\s+/,$line);
    		my $match = 0;
    		foreach my $ign(@$ignore_contigs){
    			if("$ign" eq "$contig"){
    				$match = 1;
    				last;
    			}
    		}
    		next if($match);
    		$contig_count++;
    	}
    close($FH);
  }
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

1;
