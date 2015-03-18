#!/bin/bash

rm -f suppl-realdata.tex

for i in figures/HCP*.pdf; do
    fname=$(basename $i)
    fname=${fname::-4}
    echo "\\begin{figure}[!ht]" >> suppl-realdata.tex
    echo "  \\centering" >> suppl-realdata.tex
    echo "  \\textbf{$fname}" >> suppl-realdata.tex
    echo "  \\includegraphics[width=.9\textwidth]{$i}" >> suppl-realdata.tex
    echo "\\end{figure}" >> suppl-realdata.tex
done
