#!/bin/bash
#PBS -l select=1:ncpus=20:mpiprocs=WIN:ngpus=4
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -N COMPacesREPN
#PBS -q gpu
#

module load gcc-8.5.0/amber22
unset CUDA_VISIBLE_DEVICES

cd $PBS_O_WORKDIR
mpirun -np WIN pmemd.cuda_SPFP.MPI -ng WIN -groupfile GROUPFILE -rem 3

