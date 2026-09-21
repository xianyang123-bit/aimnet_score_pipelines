#!/bin/bash
set -euo pipefail
package=/home/xianyang/aimnet2_score_pipelines
run=$package/work/L1234-fixed-receptor-eint-20260921
python=/home/xianyang/miniconda3/envs/t3-aimnet/bin/python
deps=$package/work/pocket-correction-20260921/deps
entry=$package/work/eint_pipeline/run.py
logs=$run/logs

pipeline=$(sbatch --parsable --job-name=eint-prepare-dock --partition=preempt --qos=preempt_cpu_qos --cpus-per-task=4 --mem=10G --time=12:00:00 --array=0-15%8 --output=$logs/pipeline-%A_%a.log --wrap="cd $package && PYTHONPATH=$deps $python $entry pipeline --run $run --worker \$SLURM_ARRAY_TASK_ID --workers 16")
score=$(sbatch --parsable --dependency=afterany:$pipeline --job-name=eint-score --partition=general --gres=gpu:1 --cpus-per-task=4 --mem=20G --time=06:00:00 --array=0-31%6 --output=$logs/score-%A_%a.log --wrap="cd $package && PYTHONPATH=$deps $python $entry score --run $run --worker \$SLURM_ARRAY_TASK_ID --workers 32")
aggregate=$(sbatch --parsable --dependency=afterany:$score --job-name=eint-aggregate --partition=preempt --qos=preempt_cpu_qos --cpus-per-task=2 --mem=8G --time=00:30:00 --output=$logs/aggregate-%j.log --wrap="cd $package && PYTHONPATH=$deps $python $entry aggregate --run $run")
printf '{"pipeline":"%s","score":"%s","aggregate":"%s"}\n' "$pipeline" "$score" "$aggregate" > $run/jobs.json
cat $run/jobs.json
