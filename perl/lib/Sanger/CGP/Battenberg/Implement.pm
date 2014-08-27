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
const my $ALLELE_COUNT_SCRIPT => q{alleleCounter};

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
const my $HETDATA_TAR_DIR => 'hetdata';
const my $HETDATA_TAR => 'hetdata.tar.gz';
const my $HET_DATA_OUTPUT => '%s_chr%s_heterozygousData.png';
const my $SUBCLONE_TAR => 'subclones.tar.gz';
const my $SUBCLONE_PNG	=> '%s_subclones_chr%d';
const my $SUBCLONE_DIR => 'subclones';
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





			#TODO Normal contamination?
			#TODO CN to VCF?
			#tarball of results
			#Zip each result type into separate tarball and tar the lot after.



	#my @commands = ($command, $bgzip, $tabix);

	#PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), \@commands, 0);

  #unlink $new_vcf;

  #PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);

  return 1;

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
