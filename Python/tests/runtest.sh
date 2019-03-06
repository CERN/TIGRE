#!/bin/bash
function runAlgs() {
   arr=("$@")
   for i in "${arr[@]}";
      do
          python ~/TIGRE/Python/tests/runscript1.py configuration1 "$i"
          python ~/TIGRE/Python/tests/runscript1.py configuration2 "$i"
          python ~/TIGRE/Python/tests/runscript1.py configuration3 "$i"
          python ~/TIGRE/Python/tests/runscript1.py configuration4 "$i"
          python ~/TIGRE/Python/tests/runscript1.py configuration5 "$i"
      done

}

array=("awasd_pocs" "asd_pocs")

runAlgs "${array[@]}"

python ~/TIGRE/Python/tests/testforlog.py log_summary Wed_Mar_06