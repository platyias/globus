#!/bin/bash
#$ -N cgibbs_single
#
#$ -cwd
#$ -S /bin/bash
#$ -l mem=4G,time=5::
#$ -o CGibbs_single.log
#$ -e CGibbs_single.err

begin=$(date +%s)
./CGibbs $1 $2 $3
end=$(date +%s)
tot_time=$(expr $end - $begin)
echo "O3 - It takes $tot_time seconds to complete this task..."
