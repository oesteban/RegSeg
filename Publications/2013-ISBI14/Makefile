SUBDIRS = Figures 

all: 00-main.pdf

00-main.pdf: 00-main.tex 00-main.bbl 00-main.aux $(SUBDIRS)
	pdflatex -shell-escape $< 
	grep Rerun 00-main.log && pdflatex $< 
	grep Rerun 00-main.log && pdflatex $< 
	touch 00-main.bbl
	touch $@


$(SUBDIRS):
	$(MAKE) -C $@

00-main.bbl: 00-main.aux 99-references.bib
	bibtex 00-main

00-main.aux: *.tex
	pdflatex -shell-escape 00-main.tex 

99-references.bib:
	echo "" > 99-references.bib
	cat ../Bibliography/*.bib >> 99-references.bib



cleanall: clean
	rm -f *.pdf *.acn *.glo *.ist *.toc *.spl *.dvi *.ps

clean:
	rm -f *.aux *.log *.out *.blg *.bbl
