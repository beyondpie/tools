library(rtracklayer)
# * human genome size
GRCh38p14fnm <- "/projects/ps-renlab2/szu/genome/gencode.v47.annotation.gtf.gz"
hg38gr <- rtracklayer::import(GRCh38p14fnm)
gene_hg38 <- hg38gr[hg38gr$type == "gene"]
r <- data.frame(
  chrom = as(seqnames(gene_hg38), "vector"),
  startFrom = start(gene_hg38),
  endTo = end(gene_hg38),
  name = gene_hg38$gene_id
) |> unique()
data.table::fwrite(r,
  file = "GRCh38.p14.simple.gene.annot.tsv",
  sep = "\t", col.names = F, row.names = F
)

hmgene_id2name <- data.frame(
  id = gene_hg38$gene_id,
  name = gene_hg38$gene_name
) |> unique()

data.table::fwrite(
  hmgene_id2name,
  file = "GRCh38.p14.id2name.tsv",
  sep = "\t", col.names = F, row.names = F
)
