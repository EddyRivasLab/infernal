# do the 7 non-filtered and 7 filtered searches with cmsearch in MPI mode on 32 CPUs
qsub -R y  -N mpi.trna.f    -o mpi.trna.f.out  -b y -cwd -V -j y -pe openmpi 32 'mpirun -mca btl self,tcp --prefix /usr/local/openmpi -np 32 cmsearch --mpi ../trna.cm 10Mb.fa > mpi.trna.f.cmsearch'
qsub -R y  -N mpi.trna.nf   -o mpi.trna.nf.out -b y -cwd -V -j y -pe openmpi 32 'mpirun -mca btl self,tcp --prefix /usr/local/openmpi -np 32 cmsearch --mpi --fil-no-hmm --fil-no-qdb ../trna.cm 1Mb.fa > mpi.trna.nf.cmsearch'

qsub -R y  -N mpi.5s.f      -o mpi.5s.f.out  -b y -cwd -V -j y -pe openmpi 32 'mpirun -mca btl self,tcp --prefix /usr/local/openmpi -np 32 cmsearch --mpi ../5s.cm 10Mb.fa > mpi.5s.f.cmsearch'
qsub -R y  -N mpi.5s.nf     -o mpi.5s.nf.out -b y -cwd -V -j y -pe openmpi 32 'mpirun -mca btl self,tcp --prefix /usr/local/openmpi -np 32 cmsearch --mpi --fil-no-hmm --fil-no-qdb ../5s.cm 1Mb.fa > mpi.5s.nf.cmsearch'

qsub -R y  -N mpi.lysine.f  -o mpi.lysine.f.out  -b y -cwd -V -j y -pe openmpi 32 'mpirun -mca btl self,tcp --prefix /usr/local/openmpi -np 32 cmsearch --mpi ../lysine.cm 10Mb.fa > mpi.lysine.f.cmsearch'
qsub -R y  -N mpi.lysine.nf -o mpi.lysine.nf.out -b y -cwd -V -j y -pe openmpi 32 'mpirun -mca btl self,tcp --prefix /usr/local/openmpi -np 32 cmsearch --mpi --fil-no-hmm --fil-no-qdb ../lysine.cm 1Mb.fa > mpi.lysine.nf.cmsearch'

qsub -R y  -N mpi.srp.f     -o mpi.srp.f.out  -b y -cwd -V -j y -pe openmpi 32 'mpirun -mca btl self,tcp --prefix /usr/local/openmpi -np 32 cmsearch --mpi ../srp.cm 10Mb.fa > mpi.srp.f.cmsearch'
qsub -R y  -N mpi.srp.nf    -o mpi.srp.nf.out -b y -cwd -V -j y -pe openmpi 32 'mpirun -mca btl self,tcp --prefix /usr/local/openmpi -np 32 cmsearch --mpi --fil-no-hmm --fil-no-qdb ../srp.cm 1Mb.fa > mpi.srp.nf.cmsearch'

qsub -R y  -N mpi.rnasep.f  -o mpi.rnasep.f.out  -b y -cwd -V -j y -pe openmpi 32 'mpirun -mca btl self,tcp --prefix /usr/local/openmpi -np 32 cmsearch --mpi ../rnasep.cm 10Mb.fa > mpi.rnasep.f.cmsearch'
qsub -R y  -N mpi.rnasep.nf -o mpi.rnasep.nf.out -b y -cwd -V -j y -pe openmpi 32 'mpirun -mca btl self,tcp --prefix /usr/local/openmpi -np 32 cmsearch --mpi --fil-no-hmm --fil-no-qdb ../rnasep.cm 1Mb.fa > mpi.rnasep.nf.cmsearch'

qsub -R y  -N mpi.ssu.f     -o mpi.ssu.f.out  -b y -cwd -V -j y -pe openmpi 32 'mpirun -mca btl self,tcp --prefix /usr/local/openmpi -np 32 cmsearch --mpi ../ssu.cm 10Mb.fa > mpi.ssu.f.cmsearch'
qsub -R y  -N mpi.ssu.nf    -o mpi.ssu.nf.out -b y -cwd -V -j y -pe openmpi 32 'mpirun -mca btl self,tcp --prefix /usr/local/openmpi -np 32 cmsearch --mpi --fil-no-hmm --fil-no-qdb ../ssu.cm 1Mb.fa > mpi.ssu.nf.cmsearch'

# LSU requires 8Gb machines (the '-q c14.q' part)
qsub -R y  -q c14.q -N mpi.lsu.f     -o mpi.lsu.f.out  -b y -cwd -V -j y -pe openmpi 32 'mpirun -mca btl self,tcp --prefix /usr/local/openmpi -np 32 cmsearch --mpi ../lsu.cm 10Mb.fa > mpi.lsu.f.cmsearch'
qsub -R y  -q c14.q -N mpi.lsu.nf    -o mpi.lsu.nf.out -b y -cwd -V -j y -pe openmpi 32 'mpirun -mca btl self,tcp --prefix /usr/local/openmpi -np 32 cmsearch --mpi --fil-no-hmm --fil-no-qdb ../lsu.cm 1Mb.fa > mpi.lsu.nf.cmsearch'
