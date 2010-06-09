# Makefile

all:    manuscript

manuscript:
	pdflatex main
	bibtex main
	pdflatex main	
	pdflatex main

preprint:
	pdflatex "\def\preprintoutput{1} \input{main.tex}"
	bibtex main
	pdflatex "\def\preprintoutput{1} \input{main.tex}"
	pdflatex "\def\preprintoutput{1} \input{main.tex}"

screen:
	pdflatex "\def\screenoutput{1} \input{main.tex}"
	bibtex main
	pdflatex "\def\screenoutput{1} \input{main.tex}"
	pdflatex "\def\screenoutput{1} \input{main.tex}"

clean:
	-rm *~ *.aux *.log main.dvi main.out main.toc main.bbl main.blg
