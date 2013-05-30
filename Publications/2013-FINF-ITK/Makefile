all: 00-main.pdf

00-main.pdf: 00-main.bbl 00-main.aux 
	pdflatex 00-main.tex 

00-main.bbl: 00-main.aux 99-references.bib
	bibtex 00-main

00-main.aux: 00-main.tex  01-abstract.tex 
	pdflatex 00-main.tex 
