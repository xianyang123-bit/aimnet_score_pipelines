#!/bin/bash
#SBATCH --job-name=setup-ligunity-casf
#SBATCH --partition=preempt
#SBATCH --qos=preempt_cpu_qos
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/home/xianyang/logs/setup-ligunity-casf-%j.out

set -eu

root=/data/user_data/xianyang/casf-2016/ligunity
envdir="$root/env"
ckptdir="$root/checkpoint"
mkdir -p "$root" "$ckptdir" "$root/hf-cache"

if [ ! -x "$envdir/bin/python" ]; then
  /home/xianyang/miniconda3/bin/conda create -y -p "$envdir" python=3.10 pip
fi

python="$envdir/bin/python"
pip="$envdir/bin/pip"

if ! "$python" -c 'import torch; assert torch.__version__.startswith("2.0")' 2>/dev/null; then
  "$pip" install --no-cache-dir 'numpy<2' \
    torch==2.0.0 --index-url https://download.pytorch.org/whl/cu118
fi

if ! "$python" -c 'import unicore' 2>/dev/null; then
  "$pip" install --no-cache-dir \
    'https://github.com/dptech-corp/Uni-Core/releases/download/0.0.3/unicore-0.0.1%2Bcu118torch2.0.0-cp310-cp310-linux_x86_64.whl'
fi

"$pip" install --no-cache-dir 'numpy<2' 'transformers==4.36.2' \
  iopath lmdb ml_collections scipy tensorboardX tqdm tokenizers rdkit pandas

checkpoint="$ckptdir/checkpoint_avg_41-50.pt"
if [ ! -s "$checkpoint" ]; then
  wget -c -O "$checkpoint" \
    'https://huggingface.co/fengb/LigUnity_VS/resolve/main/protein_ranking_vs/checkpoint_avg_41-50.pt?download=true'
fi

HF_HOME="$root/hf-cache" "$python" -c \
  'from transformers import AutoTokenizer,AutoModelForMaskedLM; n="facebook/esm2_t12_35M_UR50D"; AutoTokenizer.from_pretrained(n); AutoModelForMaskedLM.from_pretrained(n); print("ESM cached")'

"$python" -c 'import torch,unicore,transformers,lmdb,rdkit; print("torch",torch.__version__,"cuda",torch.version.cuda,"transformers",transformers.__version__,"rdkit",rdkit.__version__)'
ls -lh "$checkpoint"
