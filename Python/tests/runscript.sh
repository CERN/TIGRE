#!/bin/bash

control_c() {
    rm $DIR/*.npy
    exit
}

trap control_c SIGINT

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# make new test file and save what test we are on (targetdir)

python $DIR/make_output_directory.py

# make configurations

python $DIR/generate_configurations.py

# decide what tests to run
# TODO: for now this is hardcoded but it should check the gitrepo for changes and run tests accordingly


algs=("FDK" "fbp" "sirt" "ossart" "cgls" "asd_pocs" "awasd_pocs")

# RUN algorithm tests
for testnr in {1..4}; do
	for alg in "${algs[@]}"; do

			python $DIR/test_config_gen.py TestSequence.test_${testnr}_${alg}

		done
	done

# run other tests
# clean
rm $DIR/*.npy
