library(RUnit)
tbl.all <- get(load("tbl.all-11492x14.RData"))
#tbl.complexes <- get(load("tbl.complexes.RData"))
tbl.complexes <- get(load("../../incoming/tbl.complexes.RData"))

checkEquals(nrow(subset(tbl.all, protein %in% subset(tbl.complexes,complex=="RIB")$uniprot)), 224)
checkEquals(nrow(subset(tbl.all,    gene %in% subset(tbl.complexes,complex=="RIB")$gene)), 224)
rib.all.genes <- subset(tbl.all,  gene %in% subset(tbl.complexes,complex=="RIB")$gene)$gene
length(rib.all.genes) # 224
length(unique(rib.all.genes)) # 77

table(tbl.complexes$complex)   # BAF RIB SRM TRR  UT   (no TRR in the (12 nov) version
                               #  13  78  48   7 540

# UT all fractions: 223 rows, 96 proteins
checkEquals(nrow(subset(tbl.all, protein %in% subset(tbl.complexes,complex=="UT")$uniprot)), 1076)
checkEquals(nrow(subset(tbl.all,    gene %in% subset(tbl.complexes,complex=="UT")$gene)), 1077)

dim(subset(tbl.all, protein %in% subset(tbl.complexes, complex=="UT")$uniprot))  # 1076 14
length(subset(tbl.all, protein %in% subset(tbl.complexes, complex=="UT")$uniprot)$gene)  # 1076
length(unique(subset(tbl.all, protein %in% subset(tbl.complexes, complex=="UT")$uniprot)$gene))  # 491

dim(subset(tbl.complexes,complex=="UT")) # 540 3

genes <- subset(tbl.all,    gene %in% subset(tbl.complexes,complex=="UT")$gene)$uniprot
length(prots) # 1076
length(genes)


dim(subset(tbl.all, protein %in% subset(tbl.complexes, complex=="BAF")$uniprot))  # 31 14
length(subset(tbl.all, protein %in% subset(tbl.complexes, complex=="BAF")$uniprot)$gene)  # 31
length(unique(subset(tbl.all, protein %in% subset(tbl.complexes, complex=="BAF")$uniprot)$gene))  # 13

