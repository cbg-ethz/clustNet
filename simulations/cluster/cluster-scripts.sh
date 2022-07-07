#!/usr/bin/env bash

cd clustering

# k_clust n_vars n_bg n_it bgedges equal_cpt_bg
# 4 20 10 20 different TRUE

bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster_benchmark_args.R 4 20 10 50 different TRUE

bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster_benchmark_args.R 4 20 5 30 different TRUE
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster_benchmark_args.R 4 20 10 30 different TRUE
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster_benchmark_args.R 4 20 15 30 different TRUE
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster_benchmark_args.R 4 20 20 30 different TRUE

bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster_benchmark_args.R 2 20 10 30 different TRUE
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster_benchmark_args.R 4 20 10 30 different TRUE
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster_benchmark_args.R 6 20 10 30 different TRUE
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster_benchmark_args.R 8 20 10 30 different TRUE

bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster_benchmark_args.R 4 10 10 30 different TRUE
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster_benchmark_args.R 4 20 10 30 different TRUE
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster_benchmark_args.R 4 30 10 30 different TRUE
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster_benchmark_args.R 4 40 10 30 different TRUE

bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster_benchmark_args.R 4 20 10 30 different FALSE
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster_benchmark_args.R 4 20 10 30 same TRUE
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < cluster_benchmark_args.R 4 20 10 30 same FALSE

