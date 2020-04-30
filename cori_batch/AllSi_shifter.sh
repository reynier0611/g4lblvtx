#!/bin/bash

#dir="/global/project/projectdirs/alice/ftorales" #replace per user, often with ~
#dir="~"
dir="/global/homes/r/reynier"
if [[ ! -e $dir/out_AllSi ]]; then
    mkdir $dir/out_AllSi
fi

source $dir/Singularity/cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/sphenix_setup.sh -n
export MYINSTALL=$dir/Singularity/install
source $dir/Singularity/cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/setup_local.sh $MYINSTALL

cd $dir/Singularity/g4lblvtx/macros/
root -b -q "Fun4All_G4_FAIR_barrel_only.C($2, \"$dir/out_AllSi/barrel_only_$1\")"
