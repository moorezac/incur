cd Fiji.app/

srun \
  --job-name=cellpose_segment \
  --partition=gpuq \
  --gres=gpu:A30:1 \
  --ntasks=1 \
  --cpus-per-task=96 \
  --time 12:00:00 \
  --mem=448GB \
  --output output-%j.log \
  ./ImageJ-linux64 \
  --headless trackmate_cellpose_segment.py
