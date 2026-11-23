#!/bin/bash
# Regenerate the tables of tigre331_vs_astra.tex from the suite results in results/.
set -e
cd "$(dirname "$0")"
W=results/win_331.json; L=results/linux_331.json
rm -f tables/*.tex
python3 make_tables.py $W tables win_ >/dev/null
python3 make_tables.py $L tables lin_ >/dev/null
for f in tables/win_*.tex; do sed -i 's/\\caption{/\\caption{Windows, RTX 4070. /' $f; done
for f in tables/lin_*.tex; do sed -i 's/\\caption{/\\caption{Linux, Quadro RTX 6000. /; s/\\label{tab:/\\label{tab:lin-/' $f; done
python3 summary.py $W results/win_330.json tables/win_summary.tex >/dev/null
python3 summary.py $L results/linux_330.json tables/lin_summary.tex >/dev/null
sed -i 's/(RTX 4070, Windows)/(Quadro RTX 6000, Linux)/; s/\\label{tab:summary}/\\label{tab:lin-summary}/' tables/lin_summary.tex
for f in tables/*summary.tex; do grep -q adjustbox $f || sed -i 's/\\begin{tabular}/\\begin{adjustbox}{max width=\\linewidth}\n\\begin{tabular}/; s/\\end{tabular}/\\end{tabular}\n\\end{adjustbox}/' $f; done
