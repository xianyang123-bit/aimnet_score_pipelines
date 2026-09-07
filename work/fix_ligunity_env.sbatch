#!/bin/bash
#SBATCH --job-name=fix-ligunity-env
#SBATCH --partition=preempt
#SBATCH --qos=preempt_cpu_qos
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/home/xianyang/logs/fix-ligunity-env-%j.out

set -eu
envdir=/data/user_data/xianyang/casf-2016/ligunity/env
"$envdir/bin/pip" install --no-cache-dir scikit-learn
"$envdir/bin/python" -c 'import sklearn; print(sklearn.__version__)'
