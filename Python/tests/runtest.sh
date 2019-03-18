#!/bin/bash
function runAlgs() {
   arr=("$@")
   for i in "${arr[@]}";
      do
          python ./tests/configurations.py configuration1 "$i"
   done

}

array=("FDK" "fbp" "sirt" "ossart" "cgls" "asd_pocs" "awasd_pocs")

runAlgs "${array[@]}"

python ./tests/testforlog.py log_summary Wed_Mar_13
