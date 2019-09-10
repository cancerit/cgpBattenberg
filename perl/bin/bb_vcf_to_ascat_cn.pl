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

use strict;
use warnings FATAL => 'all';

use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Pod::Usage qw(pod2usage);

pod2usage(2) if (@ARGV != 1);

my $vcf = $ARGV[0];

my $count = 1;

my $z = new IO::Uncompress::Gunzip($vcf, MultiStream=>1) or die "gunzip failed: $GunzipError\n";
while(my $line = <$z>) {
  next if($line =~ m/^#/);
  chomp $line;
  my ($chr, $start_0, $info, $norm, $tum) = (split /\t/, $line)[0,1,7,9,10];
  my $start = $start_0+1;
  my ($end) = $info =~ m/END=([[:digit:]]+)/;
  my ($wt_total, $wt_minor) = $norm =~ m/^\.\/\.\:([[:digit:]]+):([[:digit:]]+)/;
  my ($mt_total, $mt_minor) = $tum =~ m/^\.\/\.\:([[:digit:]]+):([[:digit:]]+)/;
  if (defined $wt_total &&  defined $wt_minor && defined $mt_total && defined $mt_minor) {
    print join(q{,}, $count++, $chr, $start, $end, $wt_total, $wt_minor, $mt_total, $mt_minor),"\n";
  }
}
close($z);

__END__

=head1 NAME

bb_vcf_to_ascat_cn.pl - Create tumour and normal copy number files used by caveman

=head1 SYNOPSIS

bb_vcf_to_ascat_cn.pl [battenberg_cn.vcf.gz]

=head1 PARAMETERS

=over 8

=item B<battenberg_cn.vcf.gz>
Compressed vcf file created during the battenberg_finalise step

=back

=head1 DESCRIPTION

B<bb_vcf_to_ascat_cn.pl> Create tumour and normal copy number files used by caveman from a gzipped battenberg copynumber vcf file

=cut



