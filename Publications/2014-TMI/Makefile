#
# Makefile for TMI2014
#
# Oscar Esteban - code@oscaresteban.es
#
# based on: http://tex.stackexchange.com/questions/122241/makefile-that-prevents-tex-directory-from-being-cluttered-with-build-output

ROOT_DIR=$(shell pwd)
BUILD_DIR=$(ROOT_DIR)/.build




all: 00-main.pdf

00-main.pdf: 00-main.bbl 00-main.aux 
	pdflatex 00-main.tex 

00-main.bbl: 00-main.aux ACWE.bib
	bibtex 00-main

00-main.aux: 00-main.tex
	mkdir -p $(BUILD_DIR)
	pdflatex -shell-escape 00-main.tex

$(BUILD_DIR)/TMI2014.tex: 00-main.tex

ACWE.bib:
	echo "" > ACWE.bib
	cat ../Bibliography/*.bib >> ACWE.bib

cleanall: clean
	rm -f *.pdf *.acn *.glo *.ist *.toc *.spl *.dvi *.ps *.ttt *.fff *.lot

clean:
	rm -f *.aux *.log *.out *.blg *.bbl