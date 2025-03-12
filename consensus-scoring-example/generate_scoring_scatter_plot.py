import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.colors as px_colors

# Load the data
df = pd.read_csv('traceback.tsv',sep='\t')
df_voroscoring = pd.read_csv('voroscoring.tsv',sep='\t', header=1)

# voroscoring is ranked already, put a rank column equal to the index
df_voroscoring['voro_rank'] = np.arange(1,len(df_voroscoring)+1)
df_voroscoring.columns = ["structure", "original_name", "md5", "voro_score", "voro_rank"]

# merge on original_name (for voroscoring) and 4_seletopclusts (for traceback)
df_merged = pd.merge(df, df_voroscoring, left_on='4_seletopclusts', right_on='original_name', how='inner')

# now the scatter plot between the two scores
emscoring_df = pd.read_csv('emscoring.tsv',sep='\t')

df_merged = pd.merge(df_merged, emscoring_df, left_on='1_emscoring', right_on='structure', how='inner')
# add the column cluster ID
df_merged['cluster'] = df_merged['4_seletopclusts'].str.split('_').str[1]
# make a scatter plot and color by cluster
colors = px_colors.qualitative.Dark24
cl_limit = 6
unique_clusters = df_merged['cluster'].unique()
for i, cluster in enumerate(unique_clusters):
    df_cluster = df_merged[df_merged['cluster'] == cluster]
    i_cluster = int(cluster)
    if i_cluster < cl_limit:
        plt.scatter(df_cluster['score'], -df_cluster['voro_score'], c=colors[i], label=f"Cluster {cluster}", marker='x')
    elif i_cluster == cl_limit:
        plt.scatter(df_cluster['score'], -df_cluster['voro_score'], c="#DDDBDA", label="Other clusters", s=20)
    else:
        # put a cross for the other clusters
        plt.scatter(df_cluster['score'], -df_cluster['voro_score'], c="#DDDBDA", s=20)

# order the legend alphabetically
#plt.legend(fontsize=12, loc='upper right')
handles, labels = plt.gca().get_legend_handles_labels()
order = [0,2,4,1,5,3]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])

plt.xlabel('HADDOCK EM score',fontsize=16,)
plt.ylabel('Voro Jury score',fontsize=16,)
# set xticks to -150 to -50
plt.xticks(np.arange(-150,-30,20))
plt.yticks(np.arange(0,0.45,0.05))
plt.savefig('scoring_scatter.png', dpi=400)
