#!/bin/bash
#$ -N cgibbs_updated
#
#$ -cwd
#$ -S /bin/bash
#$ -l mem=1G,time=5::
#$ -o CGibbs.log
#$ -e CGibbs.err
#$ -t 12-12


qsub run_one.sh $SGE_TASK_ID /ifs/scratch/c2b2/dv_lab/js5084/GLOBUS_DB/PNNL/CGibbs/metafiles/metafile_updated_3.txt 2
