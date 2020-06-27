## compare-GSNAP-kallisto-counts.R by Rohan Maddamsetti.

## do correlation analysis between GSNAP raw counts and kallisto est_counts.

library(tidyverse)

## These are the abundances for each strain, using itself as reference genome.
input.files <- data.frame(file=list.files("../results/GSNAP-quantification-comparison")) %>%
    filter(str_detect(.$file,"abundance.tsv"))

make.nice.kallisto.helper <- function(file.df) {
    # input is a one row dataframe, which is the name of the abundance file.
    filename <- file.df$file
    my.strain <- str_replace(filename,"[:alpha:]{3,4}-abundance.tsv","")
    my.sample <- str_replace(filename,"-abundance.tsv","")
    my.df <- read.csv(file.path("../results/GSNAP-quantification-comparison",filename),
                      sep='\t') %>%
        mutate(strain = my.strain) %>%
        mutate(sample = my.sample) %>%
        rename(gene_id = target_id) %>%
        ## trim off _mRNA suffix.
        mutate(gene_id = str_replace(gene_id,"_mRNA",""))
    return(my.df)
}

## make a big dataframe of the kallisto est_counts.
kallisto.counts.df <- input.files %>%
    split(.$file) %>%
    map_dfr(make.nice.kallisto.helper)

GSNAP.raw.reads.df <- read.csv("../results/GSNAP-quantification-comparison/GSNAP_raw_counts.csv") %>%
    mutate(sample = str_c("HMY",sample))

comparison.df <- inner_join(kallisto.counts.df,GSNAP.raw.reads.df)

cor(comparison.df$est_counts,comparison.df$raw_counts,method="pearson")
##cor(test$est_counts,test$raw_counts,method="kendall")
##cor(test$est_counts,test$raw_counts,method="spearman")

comparison.plot <- ggplot(comparison.df,aes(x=raw_counts,y=est_counts)) +
    facet_wrap(.~sample) +
    geom_point() +
    theme_classic() +
    ylab("kallisto est_counts") +
    xlab("GSNAP raw read counts")

ggsave("../results/GSNAP-quantification-comparison/kallisto-GSNAP-correlations.pdf",comparison.plot,w=11)

calc.sample.correlations <- function(df) {
    data.frame(sample = unique(df$sample),
               pearson.cor = cor(df$est_counts, df$raw_counts,method="pearson"),
               spearman.cor = cor(df$est_counts, df$raw_counts,method="spearman"))
}

## calculate correlations, conditional on sample.
correlations <- comparison.df %>%
    split(.$sample) %>%
    map_dfr(calc.sample.correlations)
write.csv(correlations,file="../results/GSNAP-quantification-comparison/kallisto-GSNAP-correlations.csv")
