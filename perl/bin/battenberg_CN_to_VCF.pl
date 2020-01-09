#!/usr/bin/perl

##########LICENCE##########
# Copyright (c) 2014-2019 Genome Research Ltd.
#
# Author: David Jones <cgphelp@sanger.ac.uk>
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

BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  push (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);

use Getopt::Long;
use Pod::Usage qw(pod2usage);

use Bio::DB::HTS;
use Try::Tiny;
use PCAP::Cli;


use Sanger::CGP::Vcf;
use Sanger::CGP::Vcf::VCFCNConverter;
use Sanger::CGP::Vcf::BamUtil;

use Sanger::CGP::Vcf;
use Sanger::CGP::Vcf::Contig;
use Sanger::CGP::Vcf::Sample;
use Sanger::CGP::Vcf::VcfProcessLog;

{
  my $opts = setup();
	my $reference = $opts->{'r'};
	$reference =~ s/\.fai$//g;

  #parse samples and contigs from the bam files or command line args.
  my ($mt_samples, $mt_sam) = parse_samples($opts, 'tumour', $reference);
  my ($wt_samples, $wt_sam) = parse_samples($opts, 'normal', $reference);
  my $contigs = parse_contigs($opts, $mt_sam, $wt_sam);
  undef $mt_sam;
  undef $wt_sam;

  die "No samples found in normal bam file (or command line)." if(scalar values %$wt_samples == 0);
  die "Multiple samples found in normal bam file." if(scalar values %$wt_samples > 1);
  die "No samples found in mutant bam file (or command line)." if(scalar values %$mt_samples == 0);
  die "Multiple samples found in mutant bam file." if(scalar values %$mt_samples > 1);

  #Setup the converter with the parsed contigs
  my $record_converter = new Sanger::CGP::Vcf::VCFCNConverter(
    -contigs => [values %$contigs]
  );
  $record_converter->extended_cn(1);

  my ($input_loc,$output_loc,$IN_FH,$OUT_FH);
  try{
    #Open in and out files (OR STDIN/STDOUT depending on what's been given at command line)
    if($opts->{'i'}){
      open($IN_FH,'<',$opts->{'i'});
      $input_loc = $opts->{'i'};
    }else{
      $IN_FH = \*STDIN;
      $input_loc = 'STDIN';
    }

    if($opts->{'o'}){
      open($OUT_FH, '>', $opts->{'o'});
      $output_loc = $opts->{'o'};
    }else{
      $OUT_FH = \*STDOUT;
      $output_loc = 'STDOUT';
    }

    #Generate header.
    my @process_logs = ();
    my $wt_sample = (values(%$wt_samples))[0];
    my $mt_sample = (values(%$mt_samples))[0];
    push @process_logs, new Sanger::CGP::Vcf::VcfProcessLog(
	        -input_vcf_source => basename($0),
          -input_vcf_ver => Sanger::CGP::Vcf->VERSION,
          -input_vcf_param => $opts,
	  );
    my $source = basename($0). '_v'. Sanger::CGP::Vcf->VERSION;
    print $OUT_FH $record_converter->generate_header($wt_sample,$mt_sample,
                                   \@process_logs, basename($reference), $source)
                                          or croak("Unable to write VCF header: $!");


    #Iterate through input and create a record for each.
    my $fai = Bio::DB::HTS::Fai->load($reference);
    while(<$IN_FH>){
      my $line = $_;
      chomp($line);
      next if $line =~ /^\s*chr\tstartpos/;
      my ($chr, $start, $end, $mt_cn_maj, $mt_cn_min, $mt_frac, $mt_cn_maj_sec, $mt_cn_min_sec, $mt_frac_sec) = (split('\s+',$line))[0,1,2,7,8,9,10,11,12];

      #Fix to work round instances where the major and minor cn
      #are set to NA
      my $mt_cn_tot = '.';
      if ($mt_cn_maj ne 'NA' && $mt_cn_min ne 'NA') {
        $mt_cn_tot = $mt_cn_maj + $mt_cn_min;
      }
      $mt_cn_min = '.' if ($mt_cn_min eq 'NA');

      my $mt_cn_tot_sec = '.';
      if ($mt_cn_maj_sec ne 'NA' && $mt_cn_min_sec ne 'NA') {
        $mt_cn_tot_sec = $mt_cn_maj + $mt_cn_min;
      }

      my $extended = {  'mt_fcf' => $mt_frac eq 'NA' ? '.' : $mt_frac,
                        'mt_tcs' => $mt_cn_tot_sec eq 'NA' ? '.' : $mt_cn_tot_sec,
                        'mt_mcs' => $mt_cn_min_sec eq 'NA' ? '.' : $mt_cn_min_sec,
                        'mt_fcs' => $mt_frac_sec eq 'NA' ? '.' : $mt_frac_sec,
                      };

      my $wt_cn_tot = 2;
      my $wt_cn_min = 1;
      $start--; # all symbolic ALTs require preceeding base padding

      die("Unrecognised format in CN segments passed: '$line'.") unless(defined($chr)
                                                                      && defined($start)
                                                                      && defined($end)
                                                                      && defined($wt_cn_tot)
                                                                      && defined($wt_cn_min)
                                                                      && defined($mt_cn_tot)
                                                                      && defined($mt_cn_min));

      my $start_allele = $fai->fetch("$chr:$start-$start");
      print $OUT_FH $record_converter->generate_record($chr,$start,$end,$start_allele,$wt_cn_tot,$wt_cn_min,$mt_cn_tot,$mt_cn_min, $extended);
    }

  }catch{
    die "Error! Reading from: |$input_loc| Writing to: |$output_loc|\n$_";
  }finally{
    close($IN_FH) if($opts->{'i'});
    close($OUT_FH) if($opts->{'o'});
  }

}

sub parse_contigs {
  my ($opts, $mt_sam, $wt_sam) = @_;
  my $contigs;
  if(defined $mt_sam) {
    $contigs = Sanger::CGP::Vcf::BamUtil->parse_contigs($mt_sam->header->text.$wt_sam->header->text,
                                                        $opts->{'rs'},
                                                        $opts->{'ra'});
  }
  else {
    open my $REF, '<', $opts->{'r'};
    while(my $line = <$REF>) {
      my ($name, $length) = (split /\t/, $line)[0,1];
      my $contig = new Sanger::CGP::Vcf::Contig(
        -name => $name,
        -length => $length,
        -assembly => $opts->{'ra'},
        -species => $opts->{'rs'}
      );
      $contigs->{$name} = $contig;
    }
    close $REF;
  }
  return $contigs;
}

sub parse_samples {
  my ($opts, $type, $reference) = @_;
  my ($param_mod, $sam, $samp_ref);
  if($type eq 'tumour') {
    $param_mod = 'm';
  }
  elsif($type eq 'normal') {
    $param_mod = 'w';
  }
  if(defined $opts->{'sb'.$param_mod}) {
    $sam = Bio::DB::HTS->new(-bam => $opts->{'sb'.$param_mod}, -fasta => $reference);
    $samp_ref = Sanger::CGP::Vcf::BamUtil->parse_samples($sam->header->text,
                                                          $opts->{$param_mod.'ss'},
                                                          $opts->{$param_mod.'sq'},
                                                          $opts->{$param_mod.'sa'},
                                                          $opts->{$param_mod.'sc'},
                                                          $opts->{$param_mod.'sd'},
                                                          $opts->{$param_mod.'sp'});
  }
  elsif(defined $opts->{$param_mod.'sn'}) {
    $samp_ref->{$opts->{$param_mod.'sn'}} = new Sanger::CGP::Vcf::Sample(
        -name =>              $opts->{$param_mod.'sn'},
        -study =>             $opts->{$param_mod.'ss'},
        -platform =>          $opts->{$param_mod.'sp'},
        -seq_protocol =>      $opts->{$param_mod.'sq'},
        -accession =>         $opts->{$param_mod.'sa'},
        -accession_source =>  $opts->{$param_mod.'sc'},
        -description =>       $opts->{$param_mod.'sd'},);
  }
  return ($samp_ref, $sam);
}

sub setup{
  my %opts;
  #Store the command used to run this script.
  $opts{'cmd'} = join " ", $0, @ARGV;
  my @random_args;
  GetOptions(
  				'h|help' => \$opts{'h'},
					'm|man' => \$opts{'m'},
					'v|version' => \$opts{'v'},
					'o|out-file=s' => \$opts{'o'},
          'f|in-file=s' => \$opts{'i'},
          'msn|sample-name-mut=s' => \$opts{'msn'},
          'mss|sample-study-mut=s' => \$opts{'mss'},
          'msa|sample-accession-mut=s' => \$opts{'msa'},
          'msc|sample-accession-source-mut=s' => \$opts{'msc'},
          'msp|seq-platform-mut=s' => \$opts{'msp'},
          'msq|sample-sequencing-protocol-mut=s' =>\$opts{'msq'},
          'msd|sample-description-mut=s' => \$opts{'msd'},
          'wsn|sample-name-mut=s' => \$opts{'wsn'},
          'wss|sample-study-norm=s' => \$opts{'wss'},
          'wsa|sample-accession-norm=s' => \$opts{'wsa'},
          'wsc|sample-accession-source-norm=s' => \$opts{'wsc'},
          'wsp|seq-platform-norm=s' => \$opts{'wsp'},
          'wsq|sample-sequencing-protocol-norm=s' =>\$opts{'wsq'},
          'wsd|sample-description-norm=s' => \$opts{'wsd'},
          'sbm|sample-bam-mut=s' => \$opts{'sbm'},
          'sbw|sample-bam-norm=s' => \$opts{'sbw'},
          'rs|reference-species=s' => \$opts{'rs'},
          'ra|reference-assembly=s' => \$opts{'ra'},
          'r|reference=s' => \$opts{'r'},
          '<>' => sub{push(@random_args,shift(@_));}
  ) or pod2usage(2);

  if(defined $opts{'v'}){
    print "Version: ".Sanger::CGP::Vcf->VERSION."\n";
    exit;
  }

  pod2usage(-verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});


  if($opts{'i'}){
    # can come from STDIN if not defined
    PCAP::Cli::file_for_reading('i', $opts{'i'});
  }
  unless(defined $opts{'msn'}) {
    PCAP::Cli::file_for_reading('sbm', $opts{'sbm'});
  }
  unless(defined $opts{'wsn'}) {
    PCAP::Cli::file_for_reading('sbw', $opts{'sbw'});
  }
  PCAP::Cli::file_for_reading('r', $opts{'r'});

  # required: direct input
  pod2usage(-message  => "\nERROR: msq|sample-sequencing-protocol-mut must be defined.\n", -verbose => 1,  -output => \*STDERR) if(exists $opts{'msq'} && ! defined $opts{'msq'});
  pod2usage(-message  => "\nERROR: wsq|sample-sequencing-protocol-norm must be defined.\n", -verbose => 1,  -output => \*STDERR) if(exists $opts{'wsq'} && ! defined $opts{'wsq'});

  return \%opts;
}


__END__

=head1 NAME

battenberg_CN_to_VCF.pl - Take CGP copy number algorithm segment output and convert to VCF format.

=head1 SYNOPSIS

battenberg_CN_to_VCF.pl [options]

    Required parameters:
      -sample-bam-mut                  -sbm  Mutant sample bam file.
      -sample-bam-norm                 -sbw  Normal sample bam file.
      -reference                       -r    Reference file
      -sample-sequencing-protocol-mut  -msq  Sample Sequencing Protocol.
      -sample-sequencing-protocol-norm -wsq  Sample Sequencing Protocol.

    Optional parameters:
      -in-file                         -i    Input file. [STDIN]
      -out-file                        -o    Output file [STDOUT].
      -reference-species               -rs   Reference species [BAM HEADER].
      -reference-assembly              -ra   Reference assembly [BAM HEADER].
      -sample-study-mut                -mss  Mut sample study.
      -sample-accession-mut            -msa  Mut sample accession [BAM HEADER].
      -sample-accession-source-mut     -msc  Mut sample accession source.
      -seq-platform-mut                -msp  Mut sequencing platform. [BAM HEADER]
      -sample-study-norm               -wss  Normal sample study.
      -sample-accession-norm           -wsa  Normal sample accession [BAM HEADER].
      -sample-accession-source-norm    -wsc  Normal sample accession source.
      -seq-platform-norm               -wsp  Normal sequencing platform [BAM HEADER].

    Other:
     -help     -h   Brief help message.
     -man      -m   Full documentation.
     -version  -v   Version information.


=head1 PARAMETERS

=over 8

=item B<-out-file>

File to write VCF format output to. If not provided STDOUT is used.

=item B<-in-file>

File to read segment information from. If not provided STDIN is used.

=item B<-sample-study-mut>

Study of mutant sample

=item B<-sample-accession-mut>

Accession of mutant sample

=item B<-sample-accession-source-mut>

Source of mutant sample accession

=item B<-sample-platform-mut>

Platform used to sequence mutant sample

=item B<-sample-sequencing-protocol-mut>

Protocol used in sequencing mutant sample

=item B<-sample-description-mut>

Description of mutant sample

=item B<-sample-study-norm>

Study of normal sample

=item B<-sample-accession-norm>

Accession of normal sample

=item B<-sample-accession-source-norm>

Source of normal sample accession

=item B<-sample-platform-norm>

Platform used to sequence normal sample

=item B<-sample-sequencing-protocol-norm>

Protocol used in sequencing normal sample

=item B<-sample-description-norm>

Description of normal sample

=item B<-sample-bam-mut>

Mutant sample bam file

=item B<-sample-bam-norm>

Normal sample bam file

=item B<-reference-species>

Reference assembly

=item B<-reference-assembly>

Reference species

=item B<-reference>

Reference file

=item B<-version>

Print version information.

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<CN_to_VCF.pl> will attempt convert the input segmented file into VCF copy number format.

=cut

