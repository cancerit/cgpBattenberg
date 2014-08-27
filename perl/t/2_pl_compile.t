# this is a catch all to ensure all modules do compile
# added as lots of 'use' functionality is dynamic in pipeline
# and need to be sure that all modules compile.
# simple 'perl -c' is unlikely to work on head scripts any more.

use strict;
use Data::Dumper;
use Test::More;
use List::Util qw(first);
use Try::Tiny qw(try catch);
use autodie qw(:all);
use File::Find;

use FindBin qw($Bin);
my $script_path = "$Bin/../bin";

use constant COMPILE_SKIP => qw();

my $perl = $^X;

my @scripts;
find({ wanted => \&build_path_set, no_chdir => 1 }, $script_path);

for(@scripts) {
  my $script = $_;
  if( first {$script =~ m/$_$/} COMPILE_SKIP ) {
    note("SKIPPING: Script with known issues: $script");
    next;
  }
  my $message = "Compilation check: $script";
  my $command = "$perl -c $script";
  my ($pid, $process);
  try {
    $pid = open $process, $command.' 2>&1 |';
    while(<$process>){};
    close $process;
    pass($message);
  }
  catch {
    fail($message);
  };
}

done_testing();

sub build_path_set {
  push @scripts, $_ if($_ =~ m/\.pl$/);
}
