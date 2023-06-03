#!/bin/bash
#-------------------------------------------------------------------------------
#  SBATCH CONFIG
#-------------------------------------------------------------------------------
## resources
#SBATCH -A xulab
#SBATCH --partition hpc3
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=80G 
#SBATCH --time 2-00:00    
## labels and outputs
#SBATCH --job-name=scenic_withoutDask-%j.out
#SBATCH --output=scenic_c1-%j.out  # %j is the unique jobID
echo "### Starting at: $(date) ###"

source activate pyscenic
#STEP 1:

output_prefix='c1_'
output_dir='/storage/htc/joshilab/Su_Li/Spencerlab/scenic_application/results/'
output_dir=${output_dir}${output_prefix}/
mkdir ${output_dir}


f_loom_path_scenic='/storage/htc/joshilab/Su_Li/Spencerlab/scenic_application/integrated.loom'

# human tfs
#f_tfs='/storage/htc/joshilab/Su_Li/GuoquanZhang/scenic_application/hs_hgnc_curated_tfs.txt'

# mouse tfs
f_tfs='/storage/htc/joshilab/Su_Li/tools_related/scenic_AuxiliaryDatasets/allTFs_mm.txt'

#pyscenic grn ${f_loom_path_scenic} ${f_tfs} -o ${output_dir}${output_prefix}adj.csv --num_workers 2

arboreto_with_multiprocessing.py ${f_loom_path_scenic} ${f_tfs} -o ${output_dir}${output_prefix}adj.csv --num_workers 2 --method grnboost2 --seed 777

echo "step 1 finished"

# STEP 2-3:

# ranking databases
# human
#f_db_names='/storage/htc/joshilab/Su_Li/tools_related/scenic_AuxiliaryDatasets/hg19-tss-centered-5kb-7species.mc9nr.feather'

# mouse
f_db_names='/storage/htc/joshilab/Su_Li/tools_related/scenic_AuxiliaryDatasets/mm9-tss-centered-5kb-7species.mc9nr.feather'

# motif databases
# human
#f_motif_path="/storage/htc/joshilab/Su_Li/tools_related/scenic_AuxiliaryDatasets/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"

# mouse
f_motif_path="/storage/htc/joshilab/Su_Li/tools_related/scenic_AuxiliaryDatasets/motifs-v9-nr.mgi-m0.001-o0.0.tbl"
pyscenic ctx ${output_dir}${output_prefix}adj.csv \
    ${f_db_names} \
    --annotations_fname ${f_motif_path} \
    --expression_mtx_fname ${f_loom_path_scenic} \
    --output ${output_dir}${output_prefix}reg.csv \
    --mask_dropouts \
    --num_workers 20 \
    --mode "custom_multiprocessing"

echo "step 2-3 finished"



echo "### Ending at: $(date) ###"
