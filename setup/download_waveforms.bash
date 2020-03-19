#!/bin/bash

# navigate to mtbench/setup
cd $(dirname ${BASH_SOURCE[0]})


URL="https://github.com/rmodrak/mtbench.git"

mkdir ../input/
cd ../input/
git clone --branch waveforms $URL

