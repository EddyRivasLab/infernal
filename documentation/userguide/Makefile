# Makefile for producing the Infernal documentation
# CVS $Id$
###########

MANPAGES = cmalign cmbuild cmscore cmsearch 

# Preparation notes.
#   - Install latex2html (version 98-1p1)
#   - latex2html.config is customized for HMMER User Guide (beware)
#   - Install PolyglotMan (version 3.0.5)
#   - need a copy of wulogo.gif for putting in HTML directory.

# pdf:
#   Make the PDF copy of the user's guide, including converted
#   man pages.
#
pdf:
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
	   rman -f LaTeX2e ../manpages/$$file.man | ../../ssdk/rmanprocess.pl >> manpages.tex;\
	done
#        for file in $(SQUIDMANS); do\
#           rman -f LaTeX2e ../../squid/Man/$$file.man | ./rmanprocess.pl >> manpages.tex;\
#       done


clean:
	-rm -f *~ main.aux main.bbl main.blg main.dvi main.log main.out main.toc

distclean:
	-rm main.pdf
	-rm manpages.tex
	make clean
