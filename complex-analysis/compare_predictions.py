"""Perform the comparisons between haddock3 and SKEMPI data."""

import glob

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats

# Convertions between 1letter fasta and 3 letters for Amino Acids
AA_L1_CONV = {
    "A": "ALA",
    "Q": "GLN",
    "W": "TRP",
    "F": "PHE",
    "Y": "TYR",
    "C": "CYS",
    "S": "SER",
    "K": "LYS",
    "D": "ASP",
    "G": "GLY",
    "L": "LEU",
}
AA_L3_CONV = {v: k for k, v in AA_L1_CONV.items()}


def extract_alascan_csv(path):
    """Parse an alascan csv.

    Parameters
    ----------
    path : str
        Path to the anascan csv to be parsed

    Returns
    -------
    alascan_dt : dict[str, dict[str, dict[str, float]]]
        Dict holding haddock3 alascan data
    """
    alascan_dt = {}
    header = None
    with open(path, "r") as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            s_ = line.strip().split("\t")
            if not header:
                header = s_
            else:
                chainid = s_[header.index("chain")]
                resid = s_[header.index("res")]
                _ori_resname = s_[header.index("ori_resname")]
                end_resname = s_[header.index("end_resname")]
                delta_score = s_[header.index("delta_score")]
                if not chainid in alascan_dt.keys():
                    alascan_dt[chainid] = {}
                alascan_dt[chainid][resid] = {AA_L3_CONV[end_resname]: float(delta_score)}
    return alascan_dt


def extract_skempi(path, query_entry, print_alascan_parameters = False):
    """Extracts information present in the SKEMPI csv.

    Parameters
    ----------
    path : str
        Path to the SKEMPI csv data
    query_entry : str
        Name of the PDB to be analysed (PDB_Chain1_Chain2)
    print_alascan_parameters : bool, optional
        Should the function print haddock3 parameters, by default False

    Returns
    -------
    data: dict[str, dict[str, dict[str, float]]]
        Dict holding SKEMPI information
    """
    data = {}
    data2 = {}
    with open(path, "r") as fin:
        # Extract header
        header = fin.readline().strip().split(";")
        # Loop over lines
        for line in fin:
            s_ = line.strip().split(";")
            entry = s_[header.index("#Pdb")]
            # Skip if not queried entry
            if entry != query_entry:
                continue
            # Extract data
            mutation = s_[header.index("Mutation(s)_PDB")]
            _temperature = s_[header.index("Temperature")]
            _wt_affinity = s_[header.index("Affinity_wt (M)")]
            _mut_affinity = s_[header.index("Affinity_mut (M)")]
            # Limit ourselves to single point mutations
            if len(mutation.split(",")) > 1:
                continue
            # Skip if not affinity data
            if any([v == "" for v in (_wt_affinity, _mut_affinity)]):
                continue

            # Disect mutated resiudes info
            _wt_resname = mutation[0]
            mut_resname = mutation[-1]
            chainid = mutation[1]
            resid = mutation[2:-1]

            # Cast to float
            temperature = float(_temperature.replace("(assumed)", ""))
            wt_affinity = float(_wt_affinity)
            mut_affinity = float(_mut_affinity)

            # Compute delta Gs (RT * ln(Kd))
            dG_wt = (8.314/4184) * temperature * np.log(wt_affinity)
            dG_mut = (8.314/4184) * temperature * np.log(mut_affinity)
            #ddG = dG_mut - dG_wt
            ddG = dG_wt - dG_mut

            # Hold data
            if not chainid in data.keys():
                data[chainid] = {}
            if not resid in data[chainid].keys():
                data[chainid][resid] = {}
            data[chainid][resid][mut_resname] = round(ddG, 3)

            if not mut_resname in data2.keys():
                data2[mut_resname] = {}
            if not chainid in data2[mut_resname].keys():
                data2[mut_resname][chainid] = {}
            data2[mut_resname][chainid][resid] = round(ddG, 3)


    if print_alascan_parameters:
        # This part is used to generate the haddock3 configuration file
        for resname in data2.keys():
            print()
            print("[alascan]")
            print(f'scan_residue = "{AA_L1_CONV[resname]}"')
            for chainid, resids in data2[resname].items():
                print(f"resdic_{chainid} = {list([int(v) for v in resids])}")

    return data


def check_data_integrity(skempi_dt, alascan_dt):
    """Check that both dict holds the same keys.

    Parameters
    ----------
    skempi_dt : dict[str, dict[str, dict[str, float]]]
        Dict holding SKEMPI information
    alascan_dt : dict[str, dict[str, dict[str, float]]]
        Dict holding haddock3 alascan data
    """
    for chainid, res_muts in skempi_dt.items():
        assert chainid in alascan_dt.keys()
        for resid, muts_dG in res_muts.items():
            assert resid in alascan_dt[chainid].keys()
            for mut, _dG in muts_dG.items():
                assert alascan_dt[chainid][resid][mut]

    for chainid, res_muts in alascan_dt.items():
        assert chainid in skempi_dt.keys()
        for resid, muts_dG in res_muts.items():
            assert resid in skempi_dt[chainid].keys()
            for mut, _dG in muts_dG.items():
                assert skempi_dt[chainid][resid][mut]


def main(hd3_rundir, skempi_entry):
    """Run the comparisons between haddock3 alascan and SKEMPI data

    Parameters
    ----------
    hd3_rundir : str
        Name of the run_dir of haddock3
    skempi_entry : str
        Name of the entry to gather data in SKEMPI
    """
    # Extract skempi data
    skempi_dt = extract_skempi("skempi_v2.csv", skempi_entry)

    # Gather all alascan module output CSV paths
    alascan_fpaths = glob.glob(f"{hd3_rundir}/*_alascan/scan_emscoring_1.csv")
    # Set parsing variable
    alascan_dt = {}
    # Loop over files
    for alascan_fpath in alascan_fpaths:
        # Skip last alascan (full interface ALA mutation)
        if "14_alascan" in alascan_fpath:
            continue
        # Extract data
        scan_dt = extract_alascan_csv(alascan_fpath)
        # Loop over extracted data to hold them in single variable
        for chainid, chain_muts in scan_dt.items():
            if not chainid in alascan_dt.keys():
                alascan_dt[chainid] = {}
            for resid, mut_dscore in chain_muts.items():
                if not resid in alascan_dt[chainid].keys():
                    alascan_dt[chainid][resid] = {}
                for mut, dscore in mut_dscore.items():
                    alascan_dt[chainid][resid][mut] = dscore

    # Verify we have all the data
    check_data_integrity(skempi_dt, alascan_dt)
    
    # Retrieve all data
    delta_haddock_score, skempi_delta_Gs, entry_names = [], [], []
    for chainid, res_muts in alascan_dt.items():
        for resid, muts_dScore in res_muts.items():
            for mut, dScore in muts_dScore.items():
                dG = skempi_dt[chainid][resid][mut]
                delta_haddock_score.append(dScore)
                skempi_delta_Gs.append(dG)
                entry_names.append(f"Chain{chainid}_Resid{resid}_{mut}")
    
    # Write data as csv
    with open("haddock3-alascan_VS_SKEMPI.csv", "w") as fout:
        fout.write("Entry Name,ΔHADDOCK score (wt-mut),ΔΔG SKEMPI (ΔGwt-ΔGmut)\n")
        for entry_name, hd3score, dG in zip(entry_names, delta_haddock_score, skempi_delta_Gs):
            fout.write(f"{entry_name},{hd3score},{dG}\n")
    
    # Compute correlation
    correlation_matrix = np.corrcoef(delta_haddock_score, skempi_delta_Gs)
    correlation = correlation_matrix[0][1]
    print(f"Correlation: {correlation:.3f}")

    # Performs linear regression
    slope, intercept, _r, _p, _std_err = stats.linregress(delta_haddock_score, skempi_delta_Gs)
    minx = min(delta_haddock_score)
    maxx = max(delta_haddock_score)
    y_xmin = slope * minx + intercept
    y_xmax = slope * maxx + intercept


    # Re-order scores by mutation types
    by_muts = {}
    for entry_name, hd3score, dG in zip(entry_names, delta_haddock_score, skempi_delta_Gs):
        mut_res = AA_L1_CONV[entry_name[-1]]
        if not mut_res in by_muts.keys():
            by_muts[mut_res] = [[], []]
        by_muts[mut_res][0].append(hd3score)
        by_muts[mut_res][1].append(dG)

    # Draw plot
    color_ramp = mpl.colormaps["tab20"]
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    for i, mut_res in enumerate(sorted(by_muts)):
        hd3scores, ddGs = by_muts[mut_res]
        ax.scatter(hd3scores, ddGs, label=mut_res, color=color_ramp(i / 20), s=70)
    # Draw linear regression line
    ax.plot([minx, maxx], [y_xmin, y_xmax], ls="dashed", lw=1, c="gray")
    #ax.scatter(delta_haddock_score, skempi_delta_Gs)
    ax.set_ylabel("SKEMPI ($\Delta G_{wt}-\Delta G_{mut}$)", fontsize=40)
    ax.set_xlabel("$\Delta$HADDOCK score (wt-mut)", fontsize=40)
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_title(
        "Correlation between [alascan] scores and SKEMPI "
        f"experimental values: {correlation:.3f}"
        )
    ax.legend(loc="best", title="Mutated to:")
    fig.tight_layout()
    #fig.show()
    fig.savefig("haddock3-alascan_VS_SKEMPI.png", dpi=350)


if __name__ == "__main__":
    main("analyse_1brs_barnase_barstar", "1BRS_A_D")
