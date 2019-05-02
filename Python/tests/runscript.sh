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


algs=("FDK" "fbp" "sirt" "ossart" "cgls" "asd_pocs" "awasd_pocs" "fista")
tests=("configuration1.npy" "configuration2.npy" "configuration3.npy" "configuration4.npy")

# RUN algorithm tests

for i in "${algs[@]}";
	do
		for j in "${tests[@]}"
			do
				python $DIR/algorithm_test.py $j $i
			done
	done

# run other tests

# assert true or false for output of test files and publish to xml
 python $DIR/test_config1.py
 python $DIR/test_config2.py
 python $DIR/test_config3.py
 python $DIR/test_config4.py
# https://stackoverflow.com/questions/11241781/python-unittests-in-jenkins


# clean
rm $DIR/*.npy
