all: 00-main.pdf

00-main.pdf: 00-main.bbl 00-main.aux 
	pdflatex 00-main.tex 

00-main.bbl: 00-main.aux ACWE.bib
	bibtex 00-main

00-main.aux: 00-main.tex
	pdflatex -shell-escape 00-main.tex 
