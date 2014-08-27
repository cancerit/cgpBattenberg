use strict;
use Data::Dumper;
use Test::More;
use List::Util qw(first);
use File::Find;
use Cwd;
use Try::Tiny qw(try finally);
use File::Spec;

use FindBin qw($Bin);
my $lib_path = "$Bin/../lib";

# Add modules here that cannot be instantiated (should be extended and have no 'new')
# or need a set of inputs - these should be tested in own test script
use constant MODULE_SKIP => qw();


my $init_cwd = getcwd;

my @modules;
try {
  chdir($lib_path);
  find({ wanted => \&build_module_set, no_chdir => 1 }, './');
} finally {
  chdir $init_cwd;
  die "The try block died with: @_\n" if(@_);
};

for my $mod(@modules) {
  use_ok($mod) or BAIL_OUT("Unable to 'use' module $mod");
}

for my $mod(@modules) {
  ok($mod->VERSION, "Check version inheritance exists ($mod)");
  if($mod->can('new')) { # only try new on things that have new defined
    new_ok($mod) unless( first {$mod eq $_} MODULE_SKIP );
  }
}

done_testing();

sub build_module_set {
  if($_ =~ m/\.pm$/) {

    my ($dir_str,$file) = (File::Spec->splitpath( $_ ))[1,2];
    $file =~ s/\.pm$//;
    my @dirs = File::Spec->splitdir( $dir_str );
    shift @dirs;
    push @modules, (join '::', @dirs).$file;
  }
}
