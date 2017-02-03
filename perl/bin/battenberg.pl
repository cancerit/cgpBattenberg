#!/usr/bin/perl

##########LICENCE##########
# Copyright (c) 2014-2016 Genome Research Ltd.
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

use FindBin;
use lib "$FindBin::Bin/../lib";

use File::Path qw( remove_tree make_path );
use File::Spec;
use Getopt::Long;
use Pod::Usage qw(pod2usage);
use List::Util qw(first);
use Const::Fast qw(const);
use File::Copy;

use PCAP::Cli;
use Sanger::CGP::Battenberg::Implement;

const my @VALID_PROCESS => qw( allelecount baflog imputefromaf
															impute combineimpute haplotypebafs
															cleanuppostbaf plothaplotypes combinebafs
															segmentphased fitcn subclones finalise );

const my @VALID_PROTOCOLS => qw( WGS WXS RNA );
const my @VALID_GENDERS => qw(XX XY L);

const my $DEFAULT_ALLELE_COUNT_MBQ => 20;
const my $DEFAULT_PLATFORM_GAMMA=>1;
const my $DEFAULT_PHASING_GAMMA=>1;
const my $DEFAULT_SEGMENTATION_GAMMA=>10;
const my $DEFAULT_CLONALITY_DIST=>0;
const my $DEFAULT_ASCAT_DIST=>1;
const my $DEFAULT_MIN_PLOIDY=>1.6;
const my $DEFAULT_MAX_PLOIDY=>4.8;
const my $DEFAULT_MIN_RHO=>0.1;
const my $DEFAULT_MAX_RHO=>1.0;
const my $DEFAULT_MIN_GOODNESS_OF_FIT=>0.63;
const my $DEFAULT_BALANCED_THRESHOLD=>0.51;
const my $DEFAULT_PROTOCOL => 'WGS';
const my $DEFAULT_PLATFORM => 'ILLUMINA';

my %index_max = ( 'allelecount' => -1,
									'baflog' => 1,
									'imputefromaf' => -1,
									'impute' => -1,
									'combineimpute' => -1,
									'haplotypebafs' => -1,
									'cleanuppostbaf' => -1,
									'plothaplotypes' => -1,
									'combinebafs' => 1,
									'segmentphased' => 1,
									'fitcn' => 1,
									'subclones' => 1,
									'finalise' => 1,
								);

{
	my $options = setup();
	Sanger::CGP::Battenberg::Implement::prepare($options);

	my $threads = PCAP::Threaded->new($options->{'threads'});
	&PCAP::Threaded::disable_out_err if(!exists $options->{'index'} && $options->{'threads'} == 1);

	# register multi  processes
	$threads->add_function('battenberg_allelecount', \&Sanger::CGP::Battenberg::Implement::battenberg_allelecount);
	$threads->add_function('battenberg_imputefromaf', \&Sanger::CGP::Battenberg::Implement::battenberg_imputefromaf);
	$threads->add_function('battenberg_runimpute', \&Sanger::CGP::Battenberg::Implement::battenberg_runimpute);
	$threads->add_function('battenberg_combineimpute', \&Sanger::CGP::Battenberg::Implement::battenberg_combineimpute);
  $threads->add_function('battenberg_haplotypebaf', \&Sanger::CGP::Battenberg::Implement::battenberg_haplotypebaf);
	$threads->add_function('battenberg_postbafcleanup', \&Sanger::CGP::Battenberg::Implement::battenberg_postbafcleanup);
	$threads->add_function('battenberg_plothaplotypes', \&Sanger::CGP::Battenberg::Implement::battenberg_plothaplotypes);

  #Now the single processes built into the battenberg flow in order of execution.
  $threads->run(($options->{'job_count'}*2), 'battenberg_allelecount', $options) if(!exists $options->{'process'} || $options->{'process'} eq 'allelecount');

	Sanger::CGP::Battenberg::Implement::battenberg_runbaflog($options) if(!exists $options->{'process'} || $options->{'process'} eq 'baflog');

	$threads->run($options->{'job_count'}, 'battenberg_imputefromaf', $options) if(!exists $options->{'process'} || $options->{'process'} eq 'imputefromaf');

	$threads->run($options->{'job_count'}, 'battenberg_runimpute', $options) if(!exists $options->{'process'} || $options->{'process'} eq 'impute');

	$threads->run($options->{'job_count'}, 'battenberg_combineimpute', $options) if(!exists $options->{'process'} || $options->{'process'} eq 'combineimpute');

	$threads->run($options->{'job_count'}, 'battenberg_haplotypebaf', $options) if(!exists $options->{'process'} || $options->{'process'} eq 'haplotypebafs');

	$threads->run($options->{'job_count'}, 'battenberg_postbafcleanup', $options) if(!exists $options->{'process'} || $options->{'process'} eq 'cleanuppostbaf');

	$threads->run($options->{'job_count'}, 'battenberg_plothaplotypes', $options) if(!exists $options->{'process'} || $options->{'process'} eq 'plothaplotypes');

	Sanger::CGP::Battenberg::Implement::battenberg_combinebafs($options) if(!exists $options->{'process'} || $options->{'process'} eq 'combinebafs');

	Sanger::CGP::Battenberg::Implement::battenberg_segmentphased($options) if(!exists $options->{'process'} || $options->{'process'} eq 'segmentphased');

	Sanger::CGP::Battenberg::Implement::battenberg_fitcopyno($options) if(!exists $options->{'process'} || $options->{'process'} eq 'fitcn');

	Sanger::CGP::Battenberg::Implement::battenberg_callsubclones($options) if(!exists $options->{'process'} || $options->{'process'} eq 'subclones');

	if(!exists $options->{'process'} || $options->{'process'} eq 'finalise'){
		Sanger::CGP::Battenberg::Implement::battenberg_finalise($options);
		Sanger::CGP::Battenberg::Implement::battenberg_cleanup($options) unless($options->{'noclean'} == 1);
	}
}

sub setup {
  my %opts;
  GetOptions(
  				'h|help' => \$opts{'h'},
					'm|man' => \$opts{'m'},
					'o|outdir=s' => \$opts{'outdir'},
					'tb|tumbam=s' => \$opts{'tumbam'},
					'nb|normbam=s' => \$opts{'normbam'},
					't|threads=i' => \$opts{'threads'},
					'i|index=i' => \$opts{'index'},
					'p|process=s' => \$opts{'process'},
					'u|thousand-genomes-loc=s' => \$opts{'1kgenloc'},
					'r|reference=s' => \$opts{'reference'},
					'e|impute-info=s' => \$opts{'impute_info'},
					'c|prob-loci=s' => \$opts{'prob_loci'},
					'g|logs=s' => \$opts{'lgs'},
					'ig|ignore-contigs-file=s' => \$opts{'ignore_file'},
					#The following are optional params with defaults as constants
					'q|min-bq-allcount=i' => \$opts{'mbq'},
					'sg|segmentation-gamma=i' => \$opts{'seg_gamma'},
					'pg|phasing-gamma=i' => \$opts{'phase_gamma'},
					'cd|clonality-distance=i' => \$opts{'clonality_dist'},
					'ad|ascat-distance=i' => \$opts{'ascat_dist'},
					'bt|balanced-threshold=f' => \$opts{'balanced_thresh'},
					'lg|platform-gamma=i' => \$opts{'plat_gamma'},
					'mp|min-ploidy=f' => \$opts{'min_ploidy'},
					'xp|max-ploidy=f' => \$opts{'max_ploidy'},
					'mr|min-rho=f' => \$opts{'min_rho'},
          'xr|max-rho=f' => \$opts{'max_rho'},
					'mg|min-goodness-of-fit=f' => \$opts{'min_goodness'},
          'rho|purity=f' => \$opts{'rho'},
          'psi|ploidy=f' => \$opts{'psi'},
          'ch|new_chr=s' => \$opts{'new_chr'},
          'po|new_pos=s' => \$opts{'new_pos'},
          'ic|new_min_cn=s' => \$opts{'new_min_cn'},
          'ac|new_maj_cn=s' => \$opts{'new_maj_cn'},
					'rs|species=s' => \$opts{'species'},
          'ra|assembly=s' => \$opts{'assembly'},
          'pr|protocol=s' => \$opts{'protocol'},
          'pl|platform=s' => \$opts{'platform'},
          'a|allele-counts=s' => \$opts{'allele-counts'},
          'ge|gender=s' => \$opts{'gender'},
          'gl|genderloci=s' => \$opts{'genderloci'},
          'j|jobs' => \$opts{'jobs'},
          'nc|noclean' => \$opts{'noclean'},
		) or pod2usage(2);

	pod2usage(-verbose => 0, -exitval => 0) if(defined $opts{'h'});
  pod2usage(-verbose => 2, -exitval => 0) if(defined $opts{'m'});

  # then check for no args:
  my $defined;
  for(keys %opts) { $defined++ if(defined $opts{$_}); }

	pod2usage(-msg  => "\nERROR: Options must be defined.\n", -verbose => 2,  -output => \*STDERR) unless($defined);

  if(!defined($opts{'noclean'})){
    $opts{'noclean'} = 0;
  }

  if(defined $opts{'allele-counts'}) {
    for my $bam_opt((qw(tumbam normbam))) {
      if(!exists $opts{$bam_opt} || !defined $opts{$bam_opt} || $opts{$bam_opt} =~ m/[.]bam$/) {
        pod2usage(-msg  => "\nERROR: When '-a' defined $bam_opt should be sample name.\n", -verbose => 1,  -output => \*STDERR);
      }
    }
  }
  else {
    delete $opts{'allele-counts'};
    PCAP::Cli::file_for_reading('tumbam',$opts{'tumbam'});
    PCAP::Cli::file_for_reading('normbam',$opts{'normbam'});
    #We should also check the bam indexes exist.
    PCAP::Cli::file_for_reading('tumour-bai',$opts{'tumbam'}.'.bai');
    PCAP::Cli::file_for_reading('normal-bai',$opts{'normbam'}.'.bai');
  }

  PCAP::Cli::file_for_reading('impute_info.txt',$opts{'impute_info'});


  if(defined $opts{'prob_loci'}) {
    PCAP::Cli::file_for_reading('prob_loci.txt',$opts{'prob_loci'});
  }
  else {
    $opts{'prob_loci'} = Sanger::CGP::Battenberg::Implement::get_mod_path().'/probloci.txt.gz';
  }

	@{$opts{'ignored_contigs'}} = ();
	if(exists ($opts{'ignore_file'}) && defined($opts{'ignore_file'})){
		$opts{'ignored_contigs'} = Sanger::CGP::Battenberg::Implement::read_contigs_from_file($opts{'ignore_file'});
	}


	delete $opts{'process'} unless(defined $opts{'process'});
  delete $opts{'index'} unless(defined $opts{'index'});

  if(defined $opts{'gender'}){
    pod2usage(-message => 'unknown gender value: '.$opts{'gender'}, -verbose => 1) unless(first {$_ eq $opts{'gender'}} @VALID_GENDERS);
    if($opts{'gender'} eq 'L') {
      die "ERROR: Gender of XY/XX must be supplied when 'allele-counts' defined\n" if(defined $opts{'allele-counts'});
      $opts{'gender'} = Sanger::CGP::Battenberg::Implement::determine_gender(\%opts);
    }
  } else {
    pod2usage(-message => 'gender not set', -verbose => 1);
  }

  if(exists $opts{'protocol'} && defined $opts{'protocol'}) {
    my $bad_prot = 1;
		$bad_prot = 0 if(first { $_ eq $opts{'protocol'} } @VALID_PROTOCOLS);
		pod2usage(-msg  => "\nERROR: Invalid pr|protocol '$opts{protocol}'.\n", -verbose => 1,  -output => \*STDERR) if($bad_prot);
  }

  my $no_of_jobs = Sanger::CGP::Battenberg::Implement::file_line_count_with_ignore($opts{'reference'},$opts{'ignored_contigs'});
  $opts{'job_count'} = $no_of_jobs;

	if(exists $opts{'process'}) {
    PCAP::Cli::valid_process('process', $opts{'process'}, \@VALID_PROCESS);

    my $jobs = 1;
    if(exists $opts{'process'} && defined $opts{'process'}) {
      if($opts{'process'} eq 'allelecount') {
        $jobs = $no_of_jobs*2;
      }
      elsif(first {$opts{'process'} eq $_} (qw(imputefromaf impute combineimpute haplotypebafs cleanuppostbaf plothaplotypes))) {
        $jobs = $no_of_jobs;
      }
      if(exists $opts{'jobs'} && defined $opts{'jobs'}) {
        print "Jobs to complete process '$opts{process}' = $jobs\n";
        exit;
      }
    }

    if(exists $opts{'index'}) {
      my $max = $index_max{$opts{'process'}};
      $max = $jobs if($max==-1);

      die "ERROR: based on reference and exclude option index must be between 1 and $max\n" if($opts{'index'} < 1 || $opts{'index'} > $max);
      PCAP::Cli::opt_requires_opts('index', \%opts, ['process']);

      die "No max has been defined for this process type\n" if($max == 0);

      PCAP::Cli::valid_index_by_factor('index', $opts{'index'}, $max, 1);
    }
  }
  elsif(exists $opts{'index'}) {
    die "ERROR: -index cannot be defined without -process\n";
  }

  # now safe to apply defaults
	$opts{'threads'} = 1 unless(defined $opts{'threads'});

	#Create the results directory in the output directory given.
	my $tmpdir = File::Spec->catdir($opts{'outdir'}, 'tmpBattenberg');
	make_path($tmpdir) unless(-d $tmpdir);
	$opts{'tmp'} = $tmpdir;
	my $resultsdir = File::Spec->catdir($opts{'tmp'}, 'results');
	make_path($resultsdir) unless(-d $resultsdir);
	#directory to store progress reports
	my $progress = File::Spec->catdir($opts{'tmp'}, 'progress');
  make_path($progress) unless(-d $progress);
	#Directory to store run logs.
	my $logs;
	if(defined $opts{'lgs'}){
	  $logs = $opts{'lgs'};
	}else{
    $logs = File::Spec->catdir($opts{'tmp'}, 'logs');
	}
	make_path($logs) unless(-d $logs);
	$opts{'logs'} = $logs;

	if($opts{'gender'} eq 'XY'){
		$opts{'is_male'} = 'TRUE';
	}else{
		$opts{'is_male'} = 'FALSE';
	}

  #If either psi or rho are defined, then both must be defined
  if (defined $opts{'psi'} | defined $opts{'rho'}) {
    unless (defined $opts{'rho'} & defined $opts{'psi'}) {
      die "ERROR: Please define values for purity and ploidy\n";
    }
  }

  #If any of these are defined, then all must be defined
  if (defined($opts{'new_chr'}) | defined($opts{'new_pos'}) | defined($opts{'new_min_cn'}) | defined($opts{'new_maj_cn'})) {
    unless (defined($opts{'new_chr'}) & defined($opts{'new_pos'}) & defined($opts{'new_min_cn'}) & defined($opts{'new_maj_cn'})) {
      die "ERROR: Please define values for new_chr, new_pos, new_min_cn and new_maj_cn\n";
    }
  }

  #Check we have only been supplied ploidy and purity OR new_chr, new_pos, new_min_cn, new_maj_cn
  if (defined $opts{'psi'} && defined($opts{'new_chr'})) {
    die "ERROR: Please define either ploidy and purity OR new_chr, new_pos, new_min_cn, new_maj_cn\n";
  }


	#Setup default values if they're not set at commandline

	$opts{'mbq'} = $DEFAULT_ALLELE_COUNT_MBQ if(!exists($opts{'mbq'}) || !defined($opts{'mbq'}));
	$opts{'seg_gamma'} = $DEFAULT_SEGMENTATION_GAMMA if(!exists($opts{'seg_gamma'}) || !defined($opts{'seg_gamma'}));
	$opts{'phase_gamma'} = $DEFAULT_PHASING_GAMMA if(!exists($opts{'phase_gamma'}) || !defined($opts{'phase_gamma'}));
	$opts{'clonality_dist'} = $DEFAULT_CLONALITY_DIST if(!exists($opts{'clonality_dist'}) || !defined($opts{'clonality_dist'}));
	$opts{'ascat_dist'} = $DEFAULT_ASCAT_DIST if(!exists($opts{'ascat_dist'}) || !defined($opts{'ascat_dist'}));
	$opts{'balanced_thresh'} = $DEFAULT_BALANCED_THRESHOLD if(!exists($opts{'balanced_thresh'}) || !defined($opts{'balanced_thresh'}));
	$opts{'plat_gamma'} = $DEFAULT_PLATFORM_GAMMA if(!exists($opts{'plat_gamma'}) || !defined($opts{'plat_gamma'}));
	$opts{'min_ploidy'} = $DEFAULT_MIN_PLOIDY if(!exists($opts{'min_ploidy'}) || !defined($opts{'min_ploidy'}));
	$opts{'max_ploidy'} = $DEFAULT_MAX_PLOIDY if(!exists($opts{'max_ploidy'}) || !defined($opts{'max_ploidy'}));
	$opts{'min_rho'} = $DEFAULT_MIN_RHO if(!exists($opts{'min_rho'}) || !defined($opts{'min_rho'}));
  $opts{'max_rho'} = $DEFAULT_MAX_RHO if(!exists($opts{'max_rho'}) || !defined($opts{'max_rho'}));
	$opts{'min_goodness'} = $DEFAULT_MIN_GOODNESS_OF_FIT if(!exists($opts{'min_goodness'}) || !defined($opts{'min_goodness'}));

	$opts{'protocol'} = $DEFAULT_PROTOCOL if(!exists($opts{'protocol'}) || !defined($opts{'protocol'}));
	$opts{'platform'} = $DEFAULT_PLATFORM if(!exists($opts{'platform'}) || !defined($opts{'platform'}));

	return \%opts;
}


__END__

=head1 NAME

battenberg.pl - Analyse aligned bam files via the battenberg algorithm to detect subclones and copynumber variations.

=head1 SYNOPSIS

battenberg.pl [options]

  Required parameters:
    -outdir                -o   Folder to output result to.
    -reference             -r   Path to reference genome index file *.fai
    -tumbam                -tb  Path to tumour bam file
                                 - when '-a' defined sample name
    -normbam               -nb  Path to normal bam file
                                 - when '-a' defined sample name
    -gender                -ge  Gender, XX, XY or L (see -gl)
    -impute-info           -e   Location of the impute info file
    -thousand-genomes-loc  -u   Location of the directory containing 1k genomes data
    -ignore-contigs-file   -ig  File containing contigs to ignore
                                - specifically male sex chromosome, mitochondria and non primary contigs

   Optional parameters:
    -allele-counts         -a   Provide a tar.gz containing the impute allele counts
    -prob-loci             -c   Location of prob_loci.txt file [included in release]
    -min-bq-allcount       -q   Minimum base quality to permit allele counting [20]
    -segmentation-gamma    -sg  Segmentation gamma [10]
    -phasing-gamma         -pg  Phasing gamma [1]
    -clonality-distance    -cd  Clonality distance [0]
    -ascat-distance        -ad  ASCAT distance [1]
    -balanced-threshold    -bt  Balanced threshold [0.51]
    -platform-gamma        -lg  Platform gamma [1]
    -min-ploidy            -mp  Min ploidy [1.6]
    -max-ploidy            -xp  Max ploidy [4.8]
    -min-rho               -mr  Min Rho [0.1]
    -max-rho               -xr  Max Rho [1.0]
    -min-goodness-of-fit   -mg  Min goodness of fit [0.63]
    -species               -rs  Reference species []
    -assembly              -ra  Reference assembly []
    -protocol              -pr  Sequencing protocol [WGS]
    -platform              -pl  Sequencing platfrom [ILLUMINA]
    -genderloci            -gl  List of gender loci, required when '-ge L' [share/gender/GRCh37d5_Y.loci]
                                - these are loci that will not present at all in a female sample

   Optional system related parameters:
    -threads           -t   Number of threads allowed on this machine (default 1)
    -logs              -g   Location to write logs (default is ./logs)

   Targeted processing (further detail under OPTIONS):
    -process  -p  Only process this step then exit, optionally set -index
    -index    -i  Optionally restrict '-p' to single job
    -jobs     -j  Declare with -p to determine how many jobs are needed for this step

   Other:
    -help     -h  Brief help message.
    -man      -m  Full documentation.

=head1 OPTIONS

=over 8

=item B<-outdir>

Directory to write output to.  During processing a temp folder will be generated in this area,
should the process fail B<only delete this if> you are unable to resume the process.

Final output files are:

 <tumour_sample>_copynumberprofile.png
 <tumour_sample>_hetdata.tar.gz
 <tumour_sample>_rafseg.tar.gz
 <tumour_sample>_second_nonroundedprofile.png
 <normal_sample>_allelecounts.tar.gz
 <tumour_sample>_allelecounts.tar.gz
 <tumour_sample>_distance.png
 <tumour_sample>_impute_input.tar.gz
 <tumour_sample>_nonroundedprofile.png
 <tumour_sample>_rho_and_psi.txt
 <tumour_sample>_subclones.tar.gz
 <tumour_sample>_battenberg_cn.vcf.gz
 <tumour_sample>_battenberg_cn.vcf.gz.tbi
 <tumour_sample>_Germline<tumour_sample>.png
 <tumour_sample>_impute_output.tar.gz
 <tumour_sample>_normal_contamination.txt
 <tumour_sample>_second_copynumberprofile.png
 <tumour_sample>_subclones.txt
 <tumour_sample>_hetbaf.tar.gz
 <tumour_sample>_logR_Baf_segmented.vcf.gz
 <tumour_sample>_logR_Baf_segmented.vcf.gz.tbi
 <tumour_sample>_other.tar.gz
 <tumour_sample>_second_distance.png
 <tumour_sample>_Tumor<tumour_sample>.png

=item B<-reference>

Path to genome.fa.fai file and the assumed location of its associated .fa file.

=item B<-tumbam>

Path to mapped, indexed, duplicate marked/removed tumour bam file.

=item B<-normbam>

Path to mapped, indexed, duplicate marked/removed normal bam file.

=item B<-ignore-contigs-file>

Path to ignore file containing a list of contigs to ignore

=item B<-is-male>

Flag the sample as male, otherwise female is assumed

=item B<-impute-info>

Path to the impute_info.txt file to be used by impute

=item B<-thousand-genomes-loc>

Directory containing the 1000genomes loci data

=item B<-prob-loci>

Location of the prob_loci.txt file

=item B< -min-bq-allcount>

Minimum base quality of alleles to be counted in allele counting step

=item B<-segmentation-gamma>

Optional battenberg input Segmentation gamma

=item B<-phasing-gamma>

Optional battenberg input Phasing gamma

=item B<-clonality-distance>

Optional battenberg input Clonality distance

=item B<-ascat-distance>

Optional battenberg input ASCAT distance

=item B<-balanced-threshold>

Optional battenberg input Balanced threshold

=item B<-platform-gamma>

Optional battenberg input Platform gamma

=item B<-min-ploidy>

Optional battenberg input Min ploidy

=item B<-max-ploidy>

Optional battenberg input Max ploidy

=item B<-min-rho>

Optional battenberg input Min Rho

=item B<-max-rho>

Optional battenberg input Max Rho

=item B<-min-goodness-of-fit>

Optional battenberg input Min goodness of fit

=item B<-species>

Reference species (default: Human)

=item B<-assembly>

Reference species assembly (default: 37)

=item B<-protocol>

Sequencing protocol (default: WGS)


=item B<-platform>

Sequencing platform (default HiSeq)

=item B<-threads>

Number of threads allowed on this machine (default 1)

=item B<-logs>

Location to write logs (default is ./logs)

=item B<-process>

Only process this step then exit, optionally set -index.  The order of steps is as follows:

    allelecount *
    baflog
    imputefromaf *
    impute *
    combineimpute *
    haplotypebafs *
    cleanuppostbaf *
    plothaplotypes *
    combinebafs
    segmentphased
    fitcn
    subclones
    finalise

'*' denotes that the step has parallel processing.  You can determine the total number by executing the command as:

    battenberg.pl ...... -p <PROCESS> -j

1 to this value can be used with '-i' on these processes.  For all other steps please use '-i 1'.

=item B<-index>

Optionally restrict '-p' to single job

=back

=head1 DESCRIPTION

B<battenberg.pl> will attempt to run all caveman steps automatically including collation of output files.

=cut

