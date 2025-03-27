#!/bin/bash
#SBATCH --partition={{resources.partition}}
#SBATCH --account={{resources.account}}
#SBATCH --time={{resources.runtime}}
#SBATCH --cpus-per-task={{resources.cpus-per-task}}

set -e
{exec_job}

