#!/bin/bash -e

# mtbench local path
WD=$(dirname ${BASH_SOURCE[0]})


# mtbech remote url
URL="https://github.com/rmodrak/mtbench.git"


function download() {
    remote=$1
    branch=$2
    basename=$3
    dirname=$4
    fullname="$3/$4"

    git clone --branch $branch $remote $fullname
    cat ${fullname}/part? > $fullname.tgz
    cd $basename
    tar -xzf $dirname.tgz
    }


#
# expected output from scripts/run_Silwal2016_sygine
#
download $URL WaveformsSilwal2016 $WD Silwal2016


#
# expected output from scripts/run_Alvizuri2018_sygine
#
#download $URL WaveformsAlvizuri2018 $WD Alvizuri2018
