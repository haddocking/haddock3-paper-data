Antibody-antigen example in which an antibody (`gevokizumab`) targets an antigen (`interleukin-1B`) at five possible epitopes, namely the ones found in PDB files 2kh2 4g6j 4g6m (the correct one) 7chz 7z4t. The example shows how to load multiple restraint information in HADDOCK 3.0. At the end of the docking, the results are analyzed to see if the correct epitope is found (which is the case here). NB: to generate the restraints please run

```bash
cd data
./generate_restraints.sh
```

making sure you have your `haddock3` environment activated. The script will generate ambiguous restraints for the five possible epitopes and archive them in the `ambig.tbl.tgz` file.