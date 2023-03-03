#!/bin/bash -e

# navigate to mtbench root directory
cd $(dirname ${BASH_SOURCE[0]})


URL="https://github.com/rmodrak/mtbench.git"

git clone --branch input $URL ./input

cd ./input
cat Alvizuri2018.tgz.part* > Alvizuri2018.tgz
cat Silwal2016.tgz.part* > Silwal2016.tgz
tar -xzf Alvizuri2018.tgz
tar -xzf Silwal2016.tgz

