for pdb in 2kh2 4g6j 4g6m 7chz 7z4t
do
haddock3-restraints active_passive_to_ambig gevokizumab_cdr.actpass epitope_${pdb}.actpass > ambig_${pdb}.tbl
done
tar -czf ambig.tbl.tgz ambig_*.tbl