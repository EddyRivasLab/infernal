Params: ./cmzasha --learn-mm 0 data/E_coli_NC_000913.fasta 
Build: release
Host: wingless.cs.washington.edu
----fastaFile: data/E_coli_NC_000913.fasta
----cmFile: not-a-cm-file (global)

After scanning, learned Markov model:
0-order Markov model:
order & count-dump list: ,0,2.28301e+06,2.35621e+06,2.35621e+06,2.28301e+06
conditional probs:
	A  = 0.246056
	C  = 0.253944
	G  = 0.253944
	U  = 0.246056

fracLetsThru=0  (# nucs lets thru=0/total nucs=9.27844e+06

CPU time: 6.59u 0.11s 00:00:06.70 Elapsed: 00:00:07
