all:
	(cd cfsqp;      make)
	(cd cm2hmm-1.0; make)

clean:
	-rm -f *~
	(cd cm2hmm-1.0; make clean)
	(cd cfsqp;      make clean)

distclean:
	-rm -f *~
	(cd cm2hmm-1.0; make distclean)
	(cd cfsqp;      make distclean)
	-rm Makefile
