library(TCGAbiolinks)

query <- GDCquery(
  project="TCGA-ESCA", 
  data.category="Simple Nucleotide Variation",
  data.type="Masked Somatic Mutation",
  access="open"
)

GDCdownload(query)
var_maf = read.maf("./TF_SNP.maf")
oncoplot(maf=var_maf)







