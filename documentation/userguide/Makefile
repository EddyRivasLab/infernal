# Makefile for producing the Infernal documentation
# SVN $Id$
#
# Preparation notes.
#   - Install latex2html (version 98-1p1)
#   - latex2html.config is customized for HMMER User Guide (beware)
#   - Install PolyglotMan (version 3.0.5)
#   - need a copy of wulogo.gif for putting in HTML directory.
###########

MANPAGES   = cmalign cmbuild cmcalibrate cmemit cmscore cmsearch cmstat
AIFIGURES  = stl9-63

# pdf:
#   Make the PDF copy of the user's guide, including converted
#   man pages.
#
pdf:
	make ai2pdf
	make manpages
	pdflatex main
	bibtex main
	pdflatex main
	pdflatex main

# manpages: 
#    Convert man pages to LaTeX chapter in User Guide.
#    uses PolyglotMan 3.0.5 "rman", and rmanprocess.pl script in ssdk.
manpages:
# write the manpages.tex "chapter" by conversion
	@echo "%% This chapter automatically generated. Do not edit." > manpages.tex
	@echo "\section{Manual pages}" >> manpages.tex
	for file in $(MANPAGES); do\
	   rman -f LaTeX2e ../manpages/$$file.man | ../../easel/devkit/rmanprocess.pl >> manpages.tex;\
	done

# ai2pdf: 
#      For any figures that we have only as .ai files, make symlinks
#      to .pdf files to fool pdflatex into importing them correctly.
ai2pdf:
	for fig in ${AIFIGURES}; do\
	   if test ! -e Figures/$$fig.pdf && test -e Figures/$$fig.ai; then (cd Figures; ln -s $$fig.ai $$fig.pdf); fi;\
	done

clean:
	-rm -f *~ main.aux main.bbl main.blg main.dvi main.log main.out main.toc

distclean:
	-rm -f *~ main.aux main.bbl main.blg main.dvi main.log main.out main.toc
	-rm -f main.pdf
	-rm -f manpages.tex
