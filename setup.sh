#!/bin/bash

########## LICENCE ##########
# Copyright (c) 2014-2017 Genome Research Ltd.
#
# Author: CancerIT <cgpit@sanger.ac.uk>
#
# This file is part of cgpBattenberg.
#
# AscatNGS is free software: you can redistribute it and/or modify it under
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
########## LICENCE ##########

done_message () {
    if [ $? -eq 0 ]; then
        echo " done."
        if [ "x$1" != "x" ]; then
            echo $1
        fi
    else
        echo " failed.  See setup.log file for error messages." $2
        echo "    Please check INSTALL file for items that should be installed by a package manager"
        exit 1
    fi
}

if [[ ($# -ne 1 && $# -ne 2) ]] ; then
  echo "Please provide an installation path such as /opt/pancan and optionally perl lib paths to allow, e.g."
  echo "  ./setup.sh /opt/myBundle"
  echo "OR all elements versioned:"
  echo "  ./setup.sh /opt/cgpBattenberg-X.X.X /opt/PCAP-X.X.X/lib/perl:/opt/alleleCount-X.X.X/lib/perl:/opt/cgpVcf-X.X.X/lib/perl"
  exit 0
fi

INST_PATH=$1

if [[ $# -eq 2 ]] ; then
  CGP_PERLLIBS=$2
fi

# get current directory
INIT_DIR=`pwd`

# cleanup inst_path
mkdir -p $INST_PATH/bin
cd $INST_PATH
INST_PATH=`pwd`
cd $INIT_DIR

# make sure that build is self contained
PERLROOT=$INST_PATH/lib/perl5
if [ -z ${CGP_PERLLIBS+x} ]; then
  export PERL5LIB="$PERLROOT"
else
  export PERL5LIB="$PERLROOT:$CGP_PERLLIBS"
fi

#create a location to build dependencies
SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR

# re-initialise log file
echo > $INIT_DIR/setup.log

# log information about this system
(
    echo '============== System information ===='
    set -x
    lsb_release -a
    uname -a
    sw_vers
    system_profiler
    grep MemTotal /proc/meminfo
    set +x
    echo; echo
) >>$INIT_DIR/setup.log 2>&1

PCAP=`perl -le 'eval "require $ARGV[0]" and print $ARGV[0]->VERSION' PCAP`
if [[ "x$PCAP" == "x" ]] ; then
  echo "PREREQUISITE: Please install PCAP-core (v1.12+) before proceeding:"
  echo "  https://github.com/ICGC-TCGA-PanCancer/PCAP-core/releases"
  exit 1;
else
  # need the leading 'v' on versions so comparison of those that don't have hotfix element behave correctly, i.e. (X.X, rather than X.X.X)
  GOOD_VER=`perl -Mversion -e "version->parse(q{v$PCAP}) >= version->parse(q{v}.q{1.12}) ? print qq{1\n} : print qq{0\n};"`
  if [[ $GOOD_VER -ne '1' ]]; then
    echo "PREREQUISITE: Please install PCAP-core (v1.12+) before proceeding:"
    echo "  https://github.com/ICGC-TCGA-PanCancer/PCAP-core/releases"
    exit 1;
  fi
fi

AC=`perl -le 'eval "require $ARGV[0]" and print $ARGV[0]->VERSION' Sanger::CGP::AlleleCount`
if [[ "x$AC" == "x" ]] ; then
  echo "PREREQUISITE: Please install alleleCount before proceeding:"
  echo "  https://github.com/cancerit/alleleCount/releases"
  exit 1;
fi

VCF=`perl -le 'eval "require $ARGV[0]" and print $ARGV[0]->VERSION' Sanger::CGP::Vcf`
if [[ "x$VCF" == "x" ]] ; then
  echo "PREREQUISITE: Please install cgpVcf (v1.3.1+) before proceeding:"
  echo "  https://github.com/cancerit/cgpVcf/releases"
  exit 1;
else
  # need the leading 'v' on versions so comparison of those that don't have hotfix element behave correctly, i.e. (X.X, rather than X.X.X)
  GOOD_VER=`perl -Mversion -e "version->parse(q{v$VCF}) >= version->parse(q{v}.q{1.3.1}) ? print qq{1\n} : print qq{0\n};"`
  if [[ $GOOD_VER -ne '1' ]]; then
    echo "PREREQUISITE: Please install cgpVcf (v1.3.1+) before proceeding:"
    echo "  https://github.com/cancerit/cgpVcf/releases"
    exit 1;
  fi
fi

## cpanm will have been installed by the pre-requisites:
CPANM=`which cpanm`
echo $CPANM

perlmods=( "File::ShareDir" "File::ShareDir::Install" )

set -e
for i in "${perlmods[@]}" ; do
  echo -n "Installing build prerequisite $i..."
  (
    set -x
    $INIT_DIR/perl/bin/cpanm -v --mirror http://cpan.metacpan.org -l $INST_PATH $i
    set +x
    echo; echo
  ) >>$INIT_DIR/setup.log 2>&1
  done_message "" "Failed during installation of $i."
done

cd $INIT_DIR/perl

echo -n "Installing Perl prerequisites ..."
if ! ( perl -MExtUtils::MakeMaker -e 1 >/dev/null 2>&1); then
    echo
    echo "WARNING: Your Perl installation does not seem to include a complete set of core modules.  Attempting to cope with this, but if installation fails please make sure that at least ExtUtils::MakeMaker is installed.  For most users, the best way to do this is to use your system's package manager: apt, yum, fink, homebrew, or similar."
fi
(
  set -x
  bin/cpanm -v --mirror http://cpan.metacpan.org --notest -l $INST_PATH/ --installdeps . < /dev/null
  set +x
) >>$INIT_DIR/setup.log 2>&1
done_message "" "Failed during installation of core dependencies."

echo -n "Installing cgpBattenberg ..."
(
  perl Makefile.PL INSTALL_BASE=$INST_PATH
  make
  make test
  make install
) >>$INIT_DIR/setup.log 2>&1
done_message "" "cgpBattenberg install failed."

# cleanup all junk
rm -rf $SETUP_DIR

echo
echo
echo "Please add the following to beginning of path:"
echo "  $INST_PATH/bin"
echo "Please add the following to beginning of PERL5LIB:"
echo "  $PERLROOT"
echo

exit 0
