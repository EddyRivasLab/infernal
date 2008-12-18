#MPI calibrations, 32 procs
cp trna.cm mpi.trna.cm
qsub -R y -N mpi.5s   -o mpi.5s.out   -b y -cwd -V -j y -pe openmpi 32 'mpirun -mca btl self,tcp --prefix /usr/local/openmpi -np 32 cmcalibrate --mpi mpi.trna.cm > mpi.trna.cmcalibrate'

cp 5s.cm mpi.5s.cm
qsub -R y -N mpi.5s   -o mpi.5s.out   -b y -cwd -V -j y -pe openmpi 32 'mpirun -mca btl self,tcp --prefix /usr/local/openmpi -np 32 cmcalibrate --mpi mpi.5s.cm > mpi.5s.cmcalibrate'

cp lysine.cm mpi.lysine.cm
qsub -R y -N mpi.5s   -o mpi.5s.out   -b y -cwd -V -j y -pe openmpi 32 'mpirun -mca btl self,tcp --prefix /usr/local/openmpi -np 32 cmcalibrate --mpi mpi.lysine.cm > mpi.lysine.cmcalibrate'

cp srp.cm mpi.srp.cm
qsub -R y -N mpi.5s   -o mpi.5s.out   -b y -cwd -V -j y -pe openmpi 32 'mpirun -mca btl self,tcp --prefix /usr/local/openmpi -np 32 cmcalibrate --mpi mpi.srp.cm > mpi.srp.cmcalibrate'

cp rnasep.cm mpi.rnasep.cm
qsub -R y -N mpi.5s   -o mpi.5s.out   -b y -cwd -V -j y -pe openmpi 32 'mpirun -mca btl self,tcp --prefix /usr/local/openmpi -np 32 cmcalibrate --mpi mpi.rnasep.cm > mpi.rnasep.cmcalibrate'

cp ssu.cm mpi.ssu.cm
qsub -R y -N mpi.5s   -o mpi.5s.out   -b y -cwd -V -j y -pe openmpi 32 'mpirun -mca btl self,tcp --prefix /usr/local/openmpi -np 32 cmcalibrate --mpi mpi.ssu.cm > mpi.ssu.cmcalibrate'

# LSU requires 8Gb machines (the '-q c14.q' part)
cp lsu.cm mpi.lsu.cm
qsub -q c14.q -R y -N mpi.5s   -o mpi.5s.out   -b y -cwd -V -j y -pe openmpi 32 'mpirun -mca btl self,tcp --prefix /usr/local/openmpi -np 32 cmcalibrate --mpi mpi.lsu.cm > mpi.lsu.cmcalibrate'
