#!/bin/bash
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

array=("sirt" "ossart" "FDK" "fbp")

runAlgs "${array[@]}"

python ~/TIGRE/Python/tests/testforlog.py log_summary Thu_Feb_28