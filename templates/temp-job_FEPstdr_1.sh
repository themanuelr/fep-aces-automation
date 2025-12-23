#!/bin/bash
#PBS -l select=1:ncpus=1:mpiprocs=1:ngpus=1
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -N COMPstdrREPN
#PBS -q gpu
#

cd $PBS_O_WORKDIR
module load gcc-8.5.0/amber22
unset CUDA_VISIBLE_DEVICES

mol_name=COMP
fep_type=stdr
TOP=PREFIX${mol_name}.top
TYPE=FEP${fep_type}
sNAME=PREFIX${mol_name}.$TYPE

for i in `seq -w 01 WIN`; do
        printf -v PREV "%02g" $(( 10#$i-1 ))
        printf -v CURR "%02g" $(( 10#$i ))
	NAME=${sNAME}.min
        pmemd.cuda_SPFP -O -p ${TOP} -i ${CURR}-min.in -c ${PREV}*.md.rst -r ${CURR}-${NAME}.rst -x ${CURR}-${NAME}.crd -o ${CURR}-${NAME}.out -ref ${PREV}*.md.rst -inf ${CURR}-${NAME}.mdinfo
	NAME=${sNAME}.md
	pmemd.cuda_SPFP -O -p ${TOP} -i ${CURR}-md.in -c ${CURR}*.min.rst -r ${CURR}-${NAME}.rst -x ${CURR}-${NAME}.crd -o ${CURR}-${NAME}.out -ref ${CURR}*.min.rst -inf ${CURR}-${NAME}.mdinfo
done

cd ../../aces/NREP
qsub job_FEPaces_1.sh

