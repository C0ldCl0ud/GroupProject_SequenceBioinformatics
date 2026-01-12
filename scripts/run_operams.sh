#!/bin/bash
# Wrapper to ensure OPERA-MS can find Racon, Minimap2, MUMmer from the host conda env

# Dynamically use the conda env that is currently active
if [ -n "$CONDA_PREFIX" ]; then
    export PATH="$CONDA_PREFIX/bin:$PATH"
fi

# Call OPERA-MS with all arguments passed
perl /operams/OPERA-MS.pl "$@"
