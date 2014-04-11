all: 00-main.pdf

00-main.pdf: 00-main.bbl 00-main.aux 
	pdflatex 00-main.tex 

00-main.bbl: 00-main.aux ACWE.bib
	bibtex 00-main

00-main.aux: 00-main.tex
	pdflatex -shell-escape 00-main.tex 

ACWE.bib:
	echo "" > ACWE.bib
	cat ../Bibliography/*.bib >> ACWE.bib

cleanall: clean
	rm -f *.pdf *.acn *.glo *.ist *.toc *.spl *.dvi *.ps *.ttt *.fff *.lot COMM_3778.txt COMM_3778.odp COMM_3778.bib COMM_3778.zip
	rm -rf COMM_3778/

clean:
	rm -f *.aux *.log *.out *.blg *.bbl package