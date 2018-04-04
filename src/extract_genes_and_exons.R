library(rtracklayer)
library(data.table)

###########
# GLOBALS #
###########

gff_file <- snakemake@input[["gff"]]
output_genes <- snakemake@output[["genes"]]
output_exons <- snakemake@output[["exons"]]
log_file <- snakemake@log[["log"]]

########
# MAIN # 
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

gff3 <- import.gff3(gff_file)

# subset the groups
gff_subset <- gff3[seqnames(gff3) %in% c("Group1.4", "Group3.5")]
#gff_subset <- gff3

# make a table of genes and exons
genes <- gff_subset[gff_subset$type == "mRNA"] #mRNA not gene because that's how the Amel2.0 annotation is formatted
gene_data <- data.table(as.character(seqnames(genes)),
           start(genes),
           end(genes),
           as.character(strand(genes)),
           genes$ID) #check the genes are denoted with ID- may need to manually edit gff

exons <- gff_subset[gff_subset$type == "CDS"] #CDS not gene
exon_data <- data.table(
    as.character(seqnames(exons)),
    start(exons),
    end(exons),
    as.character(strand(exons)),
    Parent = unlist(exons$ID))

# name the exons
exon_data[, rank := 1:length(V1), by = Parent]
exon_data[, Name := paste(Parent, rank, sep = ":")]
exon_data[, c("rank", "Parent") := NULL]

# write data
fwrite(gene_data, output_genes, sep = "\t", col.names = FALSE)
fwrite(exon_data, output_exons, sep = "\t", col.names = FALSE)

# log session
sessionInfo()
