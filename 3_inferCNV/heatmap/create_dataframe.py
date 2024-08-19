import numpy as np
import pandas as pd

df = pd.read_csv("./cnv_regions_allType_sum_percent.txt", header=0, index_col=0, sep="\t")
df = df.loc[df.Type=="deletion", :]
df = df.loc[:, ["Chr", "Percent"]]
grouped = df.groupby(by=df.index)
total_arr = []
for group in grouped:
	group[1].index = group[1].Chr
	sub_df = pd.DataFrame({group[0]:group[1].Percent}, index=group[1].index)
	total_arr.append(sub_df)
total_df = pd.concat(total_arr, axis=1, join="outer")
total_df = total_df.fillna(0)
my_index = [
"chr1p", "chr1q",
"chr2p", "chr2q",
"chr3p", "chr3q",
"chr4p", "chr4q",
"chr5p", "chr5q",
"chr6p", "chr6q",
"chr7p", "chr7q",
"chr8p", "chr8q",
"chr9p", "chr9q",
"chr10p", "chr10q",
"chr11p", "chr11q",
"chr12p", "chr12q",
"chr13p", "chr13q",
"chr14p", "chr14q",
"chr15p", "chr15q",
"chr16p", "chr16q",
"chr17p", "chr17q",
"chr18p", "chr18q",
"chr19p", "chr19q",
"chr20p", "chr20q",
"chr21p", "chr21q",
"chr22p", "chr22q"
]
my_columns = ["UBD+ Epithelial cells", "GAS2L3+ Epithelial cells", "GSTA1+ Epithelial cells", "KRT6C+ Epithelial cells"]
total_df = total_df.loc[my_index, my_columns]

total_df.to_csv("./heatmap_deletion.txt", header=True, index=True, sep="\t")

