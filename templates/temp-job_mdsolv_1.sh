#!/bin/bash
#PBS -l select=1:ncpus=1:mpiprocs=1:ngpus=1
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -N COMPsolv
#PBS -q gpu
#

cd $PBS_O_WORKDIR
module load gcc-8.5.0/amber22
unset CUDA_VISIBLE_DEVICES

mol_name=COMP

pmemd.cuda_SPFP -O -p ../${mol_name}.top -i 01-min.in -c ../${mol_name}.rst -r min.rst  -x min.crd  -o min.out
pmemd.cuda_SPFP -O -p ../${mol_name}.top -i 02-md.in -c min.rst     -r md-0.rst -x md-0.crd -o md-0.out
pmemd.cuda_SPFP -O -p ../${mol_name}.top -i 03-md.in   -c md-0.rst    -r md-1.rst   -x md-1.crd   -o md-1.out
pmemd.cuda_SPFP -O -p ../${mol_name}.top -i 04-md.in   -c md-1.rst    -r 00-${mol_name}.md.rst   -x md.crd   -o md.out
cp 00-${mol_name}.md.rst ../

cd ../
qsub job_FEPstdr_1.sh
cd ../rep1
qsub job_FEPstdr_1.sh
cd ../rep2
qsub job_FEPstdr_1.sh
cd ../rep3
qsub job_FEPstdr_1.sh

