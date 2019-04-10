#!/bin/bash

control_c() {
    exit
}

trap control_c SIGINT

function runAlgs() {
   arr=("$@")
   for i in "${arr[@]}";
      do
          python ~/TIGRE/Python/tests/configurations.py configuration1 "$i"
          python ~/TIGRE/Python/tests/configurations.py configuration2 "$i"
          python ~/TIGRE/Python/tests/configurations.py configuration3 "$i"
          python ~/TIGRE/Python/tests/configurations.py configuration4 "$i"
          python ~/TIGRE/Python/tests/configurations.py configuration5 "$i"
      done

}
array=("FDK" "fbp" "sirt" "ossart" "cgls" "asd_pocs" "awasd_pocs")

runAlgs "${array[@]}"

python ~/TIGRE/Python/tests/testforlog.py log_summary Wed_Mar_13

done