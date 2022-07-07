#!/usr/bin/env bash

bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < benchmark.R

bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < benchmark_kclust1.R
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < benchmark_kclust2.R
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < benchmark_kclust3.R
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < benchmark_kclust4.R

bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < benchmark_nbg1.R
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < benchmark_nbg2.R
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < benchmark_nbg3.R
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < benchmark_nbg4.R

bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < benchmark_qual-cpts1.R
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < benchmark_qual-cpts2.R
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < benchmark_qual-cpts3.R
bsub -n 100 -W 120:00 -R "rusage[mem=2048]" < benchmark_qual-cpts4.R


