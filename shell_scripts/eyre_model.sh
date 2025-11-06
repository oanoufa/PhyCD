#!/bin/bash
time="24:00:00"
mem="256G"

source /nfs/research/goldman/anoufa/shell_scripts/params.sh

# R VERSION
# Load R module only if Rscript is not available
# if ! command -v Rscript &> /dev/null; then
#   module load r/4.5.1
# fi

# export R_LIBS_USER="$HOME/R/library"
# mkdir -p "$R_LIBS_USER"

# wrap_cmd="Rscript \
#   /nfs/research/goldman/anoufa/src/test_eyre_model/eyre_model.R"
# sbatch -J eyre -t $time --mem=$mem \
#   -o /nfs/research/goldman/anoufa/data/out_err/eyre.out \
#   -e /nfs/research/goldman/anoufa/data/out_err/eyre.err \
#   --wrap="$wrap_cmd"

# PYTHON VERSION
wrap_cmd="${pyenv_path} \
  /nfs/research/goldman/anoufa/src/test_eyre_model/eyre_model.py"
sbatch -J eyre_py -t $time --mem=$mem \
  -o /nfs/research/goldman/anoufa/data/out_err/eyre_py.out \
  -e /nfs/research/goldman/anoufa/data/out_err/eyre_py.err \
  --wrap="$wrap_cmd"