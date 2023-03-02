f <- "~/github/erythro-dia/incoming/tbl.all-dia-rna.38662-14.RData"
file.exists(f)
tbl.all <- get(load(f))
dim(tbl.all)  # [1] 38662    14
tbl.all[1:5, 1:10]
wdth(40); colnames(tbl.all); wdth(90)
  #  [1] "protein"      "gene"
  #  [3] "species"      "fraction"
  #  [5] "peptideCount" "D0"
  #  [7] "D2"           "D4"
  #  [9] "D6"           "D8"
  # [11] "D10"          "D11"
  # [13] "D12"          "D14"
mtx <- as.matrix(tbl.all[, 6:14])
rownames(mtx) <- tbl.all$gene

x <- hclust(dist(mtx))
quartz(width=24, height=12)

plot(x)  # height goes 0:65.  3 branches if cut at 40
height <- 7e7
# table(cutree(x, h=height))  #    1    2    3
#                             # 1868  210   31
clusters <- cutree(x, h=height)
c.spc24 <- clusters[["SPC24"]]
length(clusters[clusters==c.spc24])
names(clusters[clusters==c.spc24])   # "WIPI2" "SPC24"

      #----------------------------
      # what do we see with WGCNA?
      #----------------------------

library(WGCNA)
enableWGCNAThreads()
mtx.adj <- adjacency(t(mtx), type="unsigned", power=5)
mtx.tom <- TOMsimilarity(mtx.adj)
dim(mtx.tom) # [1] 3600 3600

