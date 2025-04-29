Protein-glycan example derived from PDB ID 1LMQ in which the default HADDOCK 2.4 procedure (`unbound_vdw_na_default.cfg`) is compared with two new, HADDOCK 3.0-specific protocols, in which RMSD (`unbound_vdw_na_clustrmsd.cfg`) and IL-RMSD (`unbound_vdw_na_clustilrmsd.cfg`) clustering are performed after the rigid-body docking stage.

When running the following command:
```bash
python docking_run_analysis.py
```
a comparison between the default HADDOCK 2.4 protocol and the new HADDOCK 3.0 run is reported, with the analysis of the quality of the different clusters.

`ilrmsd_scatter.png` shows the relationship between HADDOCK3 score and Interface-ligand RMSD for the best 10 clusters of the `unbound_vdw_na_clustrmsd.cfg` HADDOCK3 run.