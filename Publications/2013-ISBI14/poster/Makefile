SUBDIRS = ../figures/

all: 00-main.pdf

00-main.pdf: 00-main.tex 00-main.blg 00-main.aux $(SUBDIRS)
	pdflatex -shell-escape $< 
	grep Rerun 00-main.log && pdflatex $< 
	grep Rerun 00-main.log && pdflatex $< 
	touch 00-main.blg
	touch $@


$(SUBDIRS):
	$(MAKE) -C $@

00-main.blg: 00-main.aux 99-references.bib
	biber 00-main

00-main.aux: *.tex beamerthemeBlue.sty beamerposter.sty
	pdflatex -shell-escape 00-main.tex 

99-references.bib:
	echo "" > 99-references.bib
	cat ../Bibliography/*.bib >> 99-references.bib

cleanall: clean
	rm -f *.pdf *.acn *.glo *.ist *.toc *.spl *.dvi *.ps *.nav *.xml *-blx.bib *.snm *.bcf

clean:
	rm -f *.aux *.log *.out *.blg *.bbl
