#! /bin/bash
#
#PBS -l walltime=00:25:00
#PBS -l nodes=1:ppn=28:gpus=1
#PBS -W group_list=newriver
#PBS -A p100_test
#PBS -q p100_normal_q

#PBS -j oe

cd $PBS_O_WORKDIR

module purge
module load gcc/5.2.0 cuda openmpi

sed ‘275s/.*/"massPartialAxHex3D_vRef0”,/‘ ellipticSolveSetupHex3D.c > foo.c
mv foo.c  massSolveSetupHex3D.c

make -j

for n in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}
do
./massMainHex3D ../../../meshes/cubeHexH025.msh $n
done



