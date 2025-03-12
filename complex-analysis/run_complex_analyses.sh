## HADDOCK related
# Download structure from PDB
wget https://files.rcsb.org/download/1BRS.pdb -O 1brs.pdb
# Extract chains A and B and remove water
pdb_selchain -A,D 1brs.pdb | pdb_delhetatm | pdb_delresname -HOH | pdb_tidy > 1brs_AD.pdb
# Run haddock3
haddock3 hd3-complex-analyses.cfg

## SKEMPI related
# Download the SKEMPI database
https://life.bsc.es/pid/skempi2/database/download/skempi_v2.csv

## Comparisons
python3 compare_predictions.py
