#!/usr/bin/perl

##########LICENCE##########
# Copyright (c) 2014 Genome Research Ltd.
#
# Author: David Jones <cgpit@sanger.ac.uk>
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

BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  push (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);

use File::Path qw( remove_tree make_path );
use File::Spec;
use Getopt::Long;
use File::Copy qw(move);
use Pod::Usage qw(pod2usage);
use Const::Fast qw(const);
use LWP::Simple qw(getstore is_success);
use Archive::Extract;
use Carp;

use Data::Dumper;

use Sanger::CGP::Battenberg;

const my $DEFAULT_URL => q{https://mathgen.stats.ox.ac.uk/impute/};
const my $IMPUTE_TGZ_PATTERN => q{ALL_1000G_phase1integrated_%s_impute.tgz};
const my $IMPUTE_LEGEND_TGZ_PATTERN => q{ALL_1000G_phase1integrated_%s_annotated_legends.tgz};
const my $IMPUTE_UNPACK_PATTERN => q{ALL_1000G_phase1integrated_%s_impute};
const my $IMPUTE_HAP_PATTERN => q{ALL_1000G_phase1integrated_%s_chr%s_impute.hap.gz};
const my $IMPUTE_LEGEND_PATTERN => q{ALL_1000G_phase1integrated_%s_chr%s_impute.legend.gz};
const my $GENETIC_MAP_PATTERN => q{genetic_map_chr%s_combined_b37.txt};
const my $IMPUTE_INFO_FILENAME => q{impute_info.txt};
const my $ONEKGEN_ALLELE_FILE_PATTERN => q{1000genomesAlleles2012_chr%s.txt};
const my $ONEKGEN_LOCI_FILE_PATTERN => q{1000genomesloci2012_chr%s.txt};


const my @THOUSAND_GENOMES_FILE_ORDER => qw(  1 2 3 4 5 6 7 8 9 10 11
																				      12 13 14 15 16 17 18 19
																				      20 21 22 23);

const my @THOUSAND_GENOMES_START => qw(  );

const my @THOUSAND_GENOMES_END => qw(  );

const my @THOUSAND_GENOMES_CHRS => qw( 1 2 3 4 5 6 7 8 9 10 11
																				      12 13 14 15 16 17 18 19
																				      20 21 22 X );

const my @IMPUTE_INFO_FILE_ORDER => qw( 1 2 3 4 5 6 7 8 9 10 11
																				12 13 14 15 16 17 18 19
																				20 21 22 X_PAR1 X_nonPAR X_PAR2);

const my @IMPUTE_START => qw{ 0 0 0 0 0 0 0 0 0 0 0
															0 0 0 0 0 0 0 0 20000000 0
															0 0 2600000 154000000};

const my @IMPUTE_END => qw{ 250000000 245000000 200000000 195000000 185000000
															175000000 160000000 150000000 145000000 140000000
															135000000 135000000 115000000 110000000 105000000
															90000000 80000000 80000000 65000000 65000000 50000000
															50000000 2700000 155000000 156000000 };

const my @IMPUTE_MALE => qw{ 1 1 1 1 1 1 1 1 1 1 1 1
															1 1 1 1 1 1 1 1 1 1 1 0 1};

const my @IMPUTE_CHRS => qw{ 1 2 3 4 5 6 7 8 9 10 11 12 13
														14 15 16 17 18 19 20 21 22 X X X};

const my %BASE_LOOKUP => ("A" => 1,
													"C" => 2,
													"G" => 3,
													"T" => 4);

const my $IMPUTE_INFO_LINE => "%s\t%s\t%s\t%s\t%d\t%d\t%d\n";
const my $ALLELE_LINE => "%d\t%d\t%d\n";
const my $LOCI_LINE => "%s\t%d\n";

my $DOWNLOAD_VERSION = "v3";

{
	my $opts = setup();
	#Download Both sets of files to directory and unpack.
	print "Download and unpack files\n";
	download_unpack_files($opts);
	#Get the date/V3 from data download.
	##copy hap.gz without replace
	my $impute_unp = sprintf($IMPUTE_UNPACK_PATTERN,$DOWNLOAD_VERSION);
	print "Move hap.gz files\n";
	my $hap_files = get_filenames_for_match(File::Spec->catdir($opts->{'tmp'},$impute_unp),$IMPUTE_HAP_PATTERN,[$DOWNLOAD_VERSION]);
	my $haps = move_files(File::Spec->catdir($opts->{'tmp'},$impute_unp),$opts->{'impdir'},$hap_files);
	##copy genetic_map_chr2_combined_b37.txt without replace
	print "Move genetic map files\n";
	my $gen_map_files = get_filenames_for_match(File::Spec->catdir($opts->{'tmp'},$impute_unp),$GENETIC_MAP_PATTERN,[]);
	my $gens = move_files(File::Spec->catdir($opts->{'tmp'},$impute_unp),$opts->{'impdir'},$gen_map_files);
	#Replace chr in legend files as we copy.
	print "Unpack and edit legend files\n";
	my $leg_files = get_filenames_for_match(File::Spec->catdir($opts->{'tmp'},$impute_unp),$IMPUTE_LEGEND_PATTERN,[$DOWNLOAD_VERSION]);
	my $legs = copy_files_with_convert(File::Spec->catdir($opts->{'tmp'},$impute_unp),$opts->{'impdir'},$leg_files);
	print "Create impute info file\n";
	create_impute_info_file($opts,$legs,$gens,$haps);
	print "Create 1000genomes files\n";
	create_one_k_genomes_files($opts,$legs);
	print "Cleaning up\n";
	cleanup($opts);
	print "Done\n";
}

sub cleanup{
  my ($opts) = @_;
  remove_tree($opts->{'tmp'}) or croak("Error removing tmporary directory $opts->{tmp}");
  return;
}

sub create_one_k_genomes_files{
	my ($opts,$legs) = @_;
	#create the 1k genomes data files.
	for(my $i=0; $i<scalar(@THOUSAND_GENOMES_FILE_ORDER); $i++){
		my $oneKGenAllFile = File::Spec->catfile($opts->{'onekdir'},
										sprintf($ONEKGEN_ALLELE_FILE_PATTERN,$THOUSAND_GENOMES_FILE_ORDER[$i]));
		my $oneKGenLociFile =  File::Spec->catfile($opts->{'onekdir'},
										sprintf($ONEKGEN_LOCI_FILE_PATTERN,$THOUSAND_GENOMES_FILE_ORDER[$i]));

		my ($ALLELE,$LOCI);
		open($ALLELE,'>',$oneKGenAllFile);
			open($LOCI,'>',$oneKGenLociFile);
				print $ALLELE "position\ta0\ta1\n" or croak ("Error writing header to allele file '$oneKGenAllFile'\n");
				#Open relevant legend file
				#Account for X having multiple sections in the legend files.
				#NB this will break if not human...
        if($THOUSAND_GENOMES_FILE_ORDER[$i] == 23){ #If we're on chromosome X
          my $j=$i;
          #use $j to iterate to the end of the list so we include all par and non par chr x positions
          for(my $j=$i; $j<scalar(@$legs); $j++){
            my $legend_file = $legs->[$j];
            my $READ;
            my $head = 0;
            my $has_snp = 0;
            open($READ,'<',$legend_file);
              while(<$READ>){
                my $line = $_;
                chomp($line);
                my ($id,$pos,$a0,$a1,$type,undef) = split(/\s+/,$line);
                if($head==0 && $line =~ m/id\s+position\s+a0/){
                  $has_snp = 1 if($type eq 'type');
                  next;
                }
                next if(($has_snp == 1  && $type ne 'SNP') || $a0 !~ m/^[ACGT]{1}$/ || $a1 !~ m/^[ACGT]{1}$/);
                $a0 = $BASE_LOOKUP{$a0};
                $a1 = $BASE_LOOKUP{$a1};
                print $ALLELE sprintf($ALLELE_LINE,$pos,$a0,$a1) or croak("Error writing 1000genomesAlleles line to '$oneKGenAllFile'\n");
                print $LOCI sprintf($LOCI_LINE,$IMPUTE_CHRS[$i],$pos) or croak("Error writing 1000genomesloci line to '$oneKGenLociFile'\n");
              }
            close($READ);
          }
        }else{
          #Only one legend file to worry about
          my $legend_file = $legs->[$i];
          my $READ;
          my $head = 0;
          my $has_snp = 0;
          open($READ,'<',$legend_file);
            while(<$READ>){
              my $line = $_;
              chomp($line);
              my ($id,$pos,$a0,$a1,$type,undef) = split(/\s+/,$line);
              if($head==0 && $line =~ m/id\s+position\s+a0/){
                $has_snp = 1 if($type eq 'type');
                next;
              }
              next if(($has_snp == 1  && $type ne 'SNP') || $a0 !~ m/^[ACGT]{1}$/ || $a1 !~ m/^[ACGT]{1}$/);
              $a0 = $BASE_LOOKUP{$a0};
              $a1 = $BASE_LOOKUP{$a1};
              print $ALLELE sprintf($ALLELE_LINE,$pos,$a0,$a1) or croak("Error writing 1000genomesAlleles line to '$oneKGenAllFile'\n");
              print $LOCI sprintf($LOCI_LINE,$IMPUTE_CHRS[$i],$pos) or croak("Error writing 1000genomesloci line to '$oneKGenLociFile'\n");
            }
          close($READ);
        }

			close($LOCI);
		close($ALLELE);
	}
	return;
}

sub create_impute_info_file{
	my ($opts,$legs,$gens,$haps) = @_;
	#Make impute info file
	my $impute_info_file = File::Spec->catfile($opts->{'impdir'},$IMPUTE_INFO_FILENAME);
	my $INFOUT;
	open ($INFOUT,'>',$impute_info_file);
		for(my $i=0;$i<scalar(@IMPUTE_END);$i++){
			print $INFOUT sprintf($IMPUTE_INFO_LINE,$IMPUTE_CHRS[$i],$legs->[$i],$gens->[$i],$haps->[$i],
														$IMPUTE_START[$i],$IMPUTE_END[$i],$IMPUTE_MALE[$i])
														or croak("Error trying to write line to impute info file '$impute_info_file'\n");
		}
	close($INFOUT);
	return;
}

sub copy_files_with_convert{
	my ($from,$to,$filenames) = @_;
	my $outfiles;
	foreach my $file(@$filenames){
		#Unpack files into the same tmp directory
		unpack_file(File::Spec->catfile($from,$file),$from,'gz');
		my $unpack_name = $file;
		$unpack_name =~ s/\.gz$//;
		my $unpack_from = File::Spec->catfile($from,$unpack_name);
		my $unpack_to = File::Spec->catfile($to,$unpack_name);
		#Open file to read
		my ($IN,$OUT);
		open($IN, '<', $unpack_from);
		#open file to write
		open($OUT, '>', $unpack_to);
			#Read each line
			while(<$IN>){
				my $line = $_;
				$line =~ s/^chr//;
				print $OUT $line or croak("Error writing chr converted legend line to $unpack_to: $!.");
			}
			#Find replace all chr with non chr
		#close read
		close($IN);
		#close write
		close($OUT);
		push(@$outfiles,$unpack_to);
	}
	return $outfiles;
}

sub move_files{
	my ($from,$to,$filenames) = @_;
	my $outfiles;
	foreach my $file(@$filenames){
		my $fromfile = File::Spec->catfile($from,$file);
		my $tofile = File::Spec->catfile($to,$file);
		croak("Error: Failed to move file $fromfile to $tofile :$!.") if(!move($fromfile,$tofile));
		push(@$outfiles,$tofile);
	}
	return $outfiles;
}

sub get_filenames_for_match{
	my ($dir,$pattern,$matchers) = @_;
	my @list = ();
	my $file_name;
	foreach my $sect(@IMPUTE_INFO_FILE_ORDER){
		push(@$matchers,$sect);
		my $file_name = sprintf($pattern,@$matchers);
		my $file = File::Spec->catfile($dir,$file_name);
		croak("Expected to find file '$file'.") if(! -e $file);
		push(@list,$file_name);
		pop(@$matchers);
	}
	return \@list;
}

sub download_unpack_files{
	my ($opts) = @_;
	my $impute_download = $opts->{'u'}.sprintf($IMPUTE_TGZ_PATTERN,$DOWNLOAD_VERSION);
	my $imputetgz = File::Spec->catfile($opts->{'tmp'}, sprintf($IMPUTE_TGZ_PATTERN,$DOWNLOAD_VERSION));
	#my $imputeout = File::Spec->catfile($opts->{'tmp'}, sprintf($IMPUTE_UNPACK_PATTERN,$DOWNLOAD_VERSION));
	download_file($impute_download,$imputetgz) if(! -e $imputetgz);
	unpack_file($imputetgz,$opts->{'tmp'},'tgz');
	return;
}

sub unpack_file{
	my ($tgz,$output_dir,$type) = @_;
	my $ae = Archive::Extract->new( archive => $tgz, type => $type );
	$ae->extract( to => $output_dir ) or croak("Error unpacking $tgz to $output_dir: ".$ae->error."\n");
	return;
}

sub download_file{
	my ($url,$file) = @_;
	croak("Error downloading $url to $file, was the version '$DOWNLOAD_VERSION' correct?\n") unless(is_success(getstore($url, $file)));
	return;
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
					'o|out-dir=s' => \$opts{'o'},
					'd|download-version=s' => \$opts{'d'},
          'u|url=s' => \$opts{'u'},
          '<>' => sub{push(@random_args,shift(@_));},
  ) or pod2usage(2);

  if(defined $opts{'v'}){
    print "Version: ".Sanger::CGP::Battenberg->VERSION."\n";
    exit;
  }

  pod2usage(-verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});
	pod2usage(-msg  => "\nERROR: Invalid inputs. Must provide o|out-dir.\n", -verbose => 1,  -output => \*STDERR) if(!exists $opts{'o'} || !defined $opts{'o'});

	#Ensure download and other directories exist, if not create it.
	my $tmpdir = File::Spec->catdir($opts{'o'}, 'tmp');
	make_path($tmpdir) unless(-d $tmpdir);

	my $imputedir = File::Spec->catdir($opts{'o'}, 'impute');
	make_path($imputedir) unless(-d $imputedir);

	my $onekdir = File::Spec->catdir($opts{'o'}, '1000genomesloci');
	make_path($onekdir) unless(-d $onekdir);

	$opts{'tmp'} = $tmpdir;
	$opts{'impdir'} = $imputedir;
	$opts{'onekdir'} = $onekdir;

	$DOWNLOAD_VERSION = $opts{'d'} if(exists($opts{'d'}) && defined($opts{'d'}));

	$opts{'u'} = $DEFAULT_URL if(!exists($opts{'u'}) || !defined($opts{'u'}));

  return \%opts;
}


__END__

=head1 NAME

download_generate_bberg_ref_files.pl - Download and create the required battenberg reference files.

=head1 SYNOPSIS

download_generate_bberg_ref_files.pl [options]

    Required parameters:
      -out-dir                  -o  Directory to output files to (see docs for structure)

    Optional parameters:
      -url                      -u  URL to download battenberg impute reference data from [see docs for default url]
      -download-version         -d  Download version (part of download name) default: [v3]

    Other:
      -help     -h   Brief help message.
      -man      -m   Full documentation.
      -version  -v   Version information.


=head1 PARAMETERS

=over 8

=item B<-out-dir>

Directory to write downloaded and manipulated files to.
Stores download in a sub tmp directory, then creates subdirectories
1000genomesloci and impute, alongside impute.info in the specified directory before deleting <your-dir>/tmp

=item B<-url>

URL to download impute files from. Default is [https://mathgen.stats.ox.ac.uk/impute/]
Files downloaded are:
ALL_1000G_phase1integrated_v3_annotated_legends.tgz
ALL_1000G_phase1integrated_v3_impute.tgz

=item B<-download-version>

version name in the download files. Default is v3.

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




