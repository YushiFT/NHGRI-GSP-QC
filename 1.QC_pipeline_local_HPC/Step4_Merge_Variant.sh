#!/bin/bash
#SBATCH -n 12
#SBATCH -N 1
#SBATCH -p general,serial_requeue,shared,xlin-lab,canstat-p01
#SBATCH -t 2880
#SBATCH --mem=50000
#SBATCH --mail-type=END
#SBATCH --mail-user=ytang@hsph.harvard.edu


module load centos6/0.0.1-fasrc01

module load python/3.6.0-fasrc01

python ../bin/Step4_Merge_Variant_1.py

cd ../temp

paste VariantQC_1.txt small_2.txt small_3.txt small_4.txt \
                      small_5.txt small_6.txt small_7.txt \
                      small_8.txt small_9.txt small_10.txt \
                      small_11.txt small_12.txt small_13.txt \
                      small_14.txt small_15.txt small_16.txt \
                      small_17.txt small_18.txt > temp.txt

rm -f smal*

python ../bin/Step4_Merge_Variant_2.py

rm -f temp.txt
