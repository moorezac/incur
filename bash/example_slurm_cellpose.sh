#! bin/bash
ssh slurm-login << 'EOF'
cd /vast/scratch/users/moore.z/projects/incur/
module load miniconda3
conda activate cellpose_sam
srun \
  --job-name=cellpose_sam \
  --partition=gpuq \
  --gres=gpu:A100:1 \
  --ntasks=1 \
  --cpus-per-task=96 \
  --time 24:00:00 \
  --mem=448GB \
  --mail-type=BEGIN,END \
  --mail-user=moore.z@wehi.edu.au \
  python python/cellpose_sam.py /vast/scratch/users/moore.z/projects/incu/emma/20240404a/data/test/composite \
    --channels 1 \
    --cellprob -1 \
    --flow 0.5 \
    --niter 2000 \
    --diameter 100 \
    --quiet
EOF
