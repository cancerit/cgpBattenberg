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

=head
const my $IMPUTE_OUTPUT_TAR_DIR => 'impute_out';
const my $IMPUTE_OUTPUT_TAR => 'impute_out.tar.gz';

const my $SUBCLONE_PNG	=> '%s_subclones_chr%d';

const my $IMPUT_OUTPUT => '%s_impute_output_chr%d_allHaplotypeInfo.txt';
=cut

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



=head

my $rho_psi_txt = get_local_rho_psi_file($root_path,$options->{'mn'});
	my ($rho,$psi,$dist) = extract_rho_psi_from_txt($rho_psi_txt);
	my $tumour_png = getLocalPngFile($root_path,[$options->{'mn'},$options->{'mn'}],$TUMOUR_PNG);
	my $normal_png = getLocalPngFile($root_path,[$options->{'mn'},$options->{'mn'}],$NORMAL_PNG);
	my $sunrise_png = getLocalPngFile($root_path,[$options->{'mn'}],$SUNRISE_PNG);
	my $copy_no_png = getLocalPngFile($root_path,[$options->{'mn'}],$COPY_NO_PNG);
	my $non_rounded_png = getLocalPngFile($root_path,[$options->{'mn'}],$NON_ROUNDED_PNG);





logs/
PD13371a_distance.png
PD13371a_impute_input_chrsample_g.txt
PD13371a_nonroundedprofile.png
PD13371a_second_nonroundedprofile.png
PD13371a_allChromosomes_heterozygousMutBAFs_haplotyped.txt
PD13371a_GermlinePD13371a.png
PD13371a_normalBAF.tab
PD13371a.logRsegmented.txt
PD13371a_normal_contamination.txt
PD13371a_rho_and_psi.txt
PD13371a_subclones.txt
PD13371a.BAFsegmented.txt
PD13371a_mutantBAF.tab
PD13371anormal_contamination.txt
PD13371a_second_copynumberprofile.png
PD13371a_TumorPD13371a.png
PD13371a_copynumberprofile.png
PD13371a_mutantLogR.tab
PD13371a_normalLogR.tab
PD13371a_second_distance.png


=cut



#rm -f /lustre/scratch112/sanger/kr2/pan_cancer_test_sets/PD13371_genome_benchmark/battenberg/battenberg_farm/tmpBattenberg/PD13371a_impute_output_chr*_*K.txt*

#const my $IMPUTE_OUTPUT => q{%s_impute_output_chr%d.txt};


	#TODO CN to VCF?
	#tarball of results
	#Zip each result type into separate tarball and tar the lot after.



	#my @commands = ($command, $bgzip, $tabix);

	#PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), \@commands, 0);

  #unlink $new_vcf;

  #PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);

  return 1;

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

=head

const my $ALLELE_COUNT_PARA => ' -b %s -o %s -l %s ';

const my @ASCAT_RESULT_FILES => qw(aberrationreliability%s.png ASCATprofile%s.png ASPCF%s.png Germline%s.png rawprofile%s.png sunrise%s.png Tumor%s.png CopyNumberCaveman%s.csv CopyNumber%s.txt SampleStatistics%s.csv);

sub ascat {
  my $options = shift;

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  my $tum_name = sanitised_sample_from_bam($options->{'tumour'});
  my $norm_name = sanitised_sample_from_bam($options->{'normal'});

  my $ascat_out = File::Spec->catdir($tmp, 'ascat');
  make_path($ascat_out) unless(-e $ascat_out);

  my $rdata = File::Spec->catfile($ascat_out,$tum_name.'.Rdata');

  my $tumcountfile = $tum_name . '.count';
  my $normcountfile = $norm_name . '.count';

  my $tumcount = File::Spec->catfile($ascat_out,$tumcountfile);
  my $normcount = File::Spec->catfile($ascat_out,$normcountfile);

  unlink $tumcount if -l $tumcount;
  unlink $normcount if -l $normcount;

  my $lntc = 'ln -s '.File::Spec->catfile(File::Spec->catdir($tmp, 'allele_count'),$tum_name.'.allct')." $tumcount";
  my $lnnc =  'ln -s '.File::Spec->catfile(File::Spec->catdir($tmp, 'allele_count'),$norm_name.'.allct')." $normcount";

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $lntc, 0);
  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $lnnc, 0);

  my $command = _which('Rscript');

  my $mod_path = dirname(abs_path($0)).'/../share';
  $mod_path = module_dir('Sanger::CGP::Ascat') unless(-e File::Spec->catdir($mod_path, 'ascat'));

  my $ascat_path = File::Spec->catdir($mod_path, 'ascat');
  my $ascat_exe = File::Spec->catfile($ascat_path,'runASCAT.R');

  $command .= " $ascat_exe";
  $command .= " $ascat_path";
  $command .= ' '.$options->{'snp_pos'};
  $command .= ' '.$options->{'snp_gc'};
  $command .= ' '.$options->{'tumour_name'};
  $command .= ' '.$tumcountfile;
  $command .= ' '.$options->{'normal_name'};
  $command .= ' '.$normcountfile;
  $command .= ' '.$options->{'gender'};
  $command .= ' '.$rdata;

  my $original_working_dir = getcwd;

  chdir($ascat_out);

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);

  chdir($original_working_dir) || die $!;

  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);

}

sub finalise {
  my $options = shift;

  my $outdir = $options->{'outdir'};
  my $tmp = $options->{'tmp'};

  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  my $tum_name = sanitised_sample_from_bam($options->{'tumour'});
  my $ascat_out = File::Spec->catdir($tmp, 'ascat');
  foreach my $f(@ASCAT_RESULT_FILES){
    my $file = sprintf($f, $tum_name);
    my $from = File::Spec->catfile($ascat_out,$file);
    my $to = File::Spec->catfile($options->{'outdir'},$file);
    move $from,$to;
  }

  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub getTumourTmp {

}

sub getNormalTmp {

}

sub getAscatTmp {

}

sub get_allele_count_file_path {
  my ($tmp,$sample_name) = @_;
  return File::Spec->catfile(File::Spec->catdir($tmp, $sample_name),'sample.allele_count');
}

sub sanitised_sample_from_bam {
  my $sample = (PCAP::Bam::sample_name(shift))[0];
  $sample =~ s/[^a-z0-9_-]/_/ig; # sanitise sample name
  return $sample;
}

sub prepare {
  my $options = shift;
  $options->{'tumour_name'} = (PCAP::Bam::sample_name($options->{'tumour'}))[0];
  $options->{'normal_name'} = (PCAP::Bam::sample_name($options->{'normal'}))[0];
  return 1;
}

sub _which {
  my $prog = shift;
  my $l_bin = $Bin;
  my $path = File::Spec->catfile($l_bin, $prog);
  $path = which($prog) unless(-e $path);
  return $path;
}
=cut

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
