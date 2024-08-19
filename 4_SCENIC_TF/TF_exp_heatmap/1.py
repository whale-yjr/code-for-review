import numpy as np
import pandas as pd

df_gene = pd.read_csv("./TF.txt", header=0, index_col=0, sep="\t")
df_data = pd.read_csv("./Epithelial_sub_TF.txt", header=0, index_col=0, sep="\t")
df = df_data.loc[df_gene.index, :]
my_columns = ["GSTA1+ Epithelial cells", "KRT6C+ Epithelial cells", "UBD+ Epithelial cells", "GAS2L3+ Epithelial cells", "Epithelial cells"]
df = df.loc[:, my_columns]

df.to_csv("./TF_exp_heatmap.txt", header=True, index=True, sep="\t")




