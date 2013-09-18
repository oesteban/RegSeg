all: 00-main.pdf

00-main.pdf: 00-main.bbl 00-main.aux 
	pdflatex 00-main.tex 

00-main.bbl: 00-main.aux 99-references.bib
	bibtex 00-main

00-main.aux: 00-main.tex
	pdflatex -shell-escape 00-main.tex 

99-references.bib:
	echo "" > 99-references.bib
	cat ../Bibliography/*.bib >> 99-references.bib

