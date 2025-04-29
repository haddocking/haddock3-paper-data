import pandas as pd
import plotly.colors as px_colors
import matplotlib.pyplot as plt


clust_capri = pd.read_csv("run_unbound_vdw_clustrmsd/11_caprieval/capri_ss.tsv", sep="\t")

colors = px_colors.qualitative.Dark24

unique_clusters = sorted(clust_capri['cluster_ranking'].unique())
# group by cluster_ranking
for i, cluster in enumerate(unique_clusters):
    if cluster > 10:
        continue
    else:
        df_cluster = clust_capri[(clust_capri['cluster_ranking'] == cluster)]
        # print a summary of the cluster in terms of dockq and ilrmsd
        df_cluster_top4 = df_cluster[df_cluster['model-cluster_ranking'] < 5]
        print(f"Cluster {cluster}:")
        print(f"  DockQ: {df_cluster_top4['dockq'].mean():.2f} +/- {df_cluster_top4['dockq'].std():.2f}")
        print(f"  Interface ligand RMSD: {df_cluster_top4['ilrmsd'].mean():.2f} +/- {df_cluster_top4['ilrmsd'].std():.2f}")
    plt.scatter(df_cluster['ilrmsd'], df_cluster['score'],
                c=colors[i], s=10, label=f"Cluster {cluster}", alpha=0.5)
    # calculate the average over the first 4 members
    avg_score = df_cluster['score'][:4].mean()
    avg_ilrmsd = df_cluster['ilrmsd'][:4].mean()
    std_score = df_cluster['score'][:4].std()
    std_ilrmsd = df_cluster['ilrmsd'][:4].std()
    plt.errorbar(avg_ilrmsd, avg_score,
                 xerr=std_ilrmsd, yerr=std_score,
                 fmt='o', color=colors[i],
                 capsize=2, capthick=2)
    
# now the other clusters
other_clusts = clust_capri[clust_capri['cluster_ranking'] > 10]

plt.scatter(other_clusts['ilrmsd'],
            other_clusts['score'],
            c="white", edgecolors="DarkSlateGrey",
            s=10, label="Other")
plt.legend(fontsize=10, loc='upper left', ncols = 2)
plt.xlabel('Interface Ligand RMSD [$\AA$]', fontsize=16)
plt.ylabel('HADDOCK score', fontsize=16)
plt.xlim(0, 15)
plt.tight_layout()
plt.savefig('ilrmsd_scatter.png', dpi=400)

# now let's analyse the default run
clust_capri_default = pd.read_csv("run_unbound_vdw_default/8_caprieval/capri_ss.tsv", sep="\t")
unique_clusters = sorted(clust_capri_default['cluster_ranking'].unique())
print(f"\n\nDefault run:")
# group by cluster_ranking
for i, cluster in enumerate(unique_clusters):
    if cluster > 10:
        continue
    else:
        df_cluster = clust_capri_default[(clust_capri_default['cluster_ranking'] == cluster)]
        df_cluster_top4 = df_cluster[df_cluster['model-cluster_ranking'] < 5]
        print(f"Cluster {cluster}:")
        print(f"  DockQ: {df_cluster_top4['dockq'].mean():.2f} +/- {df_cluster_top4['dockq'].std():.2f}")
        print(f"  Interface ligand RMSD: {df_cluster_top4['ilrmsd'].mean():.2f} +/- {df_cluster_top4['ilrmsd'].std():.2f}")
