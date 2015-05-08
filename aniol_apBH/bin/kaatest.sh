#!/bin/bash

if [[  ! -d ./bin  ]]; then
    cd ../
    if [[ ! -d ./bin ]]; then
        echo "Error: it must be executed from the main, or bin directory"
        exit
    fi
fi
#
#  Now let shell know where to find important executables
#
main=`pwd`
dirbin=$main/bin
pydir=$main/../pythia-pgs/src
pgsdir=$pydir
ERAdir=$main/../ExRootAnalysis
MAdir=$main/../MadAnalysis
webbin=$dirbin
td=$main/../td
web=0

echo $$ >> myprocid

echo p $1
echo n $2
echo t $3
echo $pydir $4

date
a=`awk '/^.*=.*nevents/{print $1}' Cards/run_card.dat`
b=`awk '/^.*=.*gridpack/{print $1}' Cards/run_card.dat`
echo Generating $a events
if [[ $b == ".true." ]]; then
    echo Generating GridPack $b
fi
echo GridPack $b
