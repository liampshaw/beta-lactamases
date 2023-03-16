theme_basic <- function () { 
  theme_bw(base_size=8) %+replace% 
    theme(
      axis.text=element_text(colour="black")
    ) %+replace% 
    theme(
      panel.grid=element_blank()
    )
}

# Genus table
metadata.df = read.csv('../data/metadata.csv',
                       header=T,
                       stringsAsFactors = F,
                       row.names = 1)
genus.counts = sort(table(metadata.df$TaxGenus))
low.abundance.genera = names(genus.counts[which(genus.counts<10)])
# assign 'other' for these genera
metadata.df$TaxGenus.simple = sapply(metadata.df$TaxGenus,
                                     function(x) ifelse(x %in% low.abundance.genera, "Other", x))

GENERA = c("Other", names(genus.counts[which(genus.counts>10)]))
REGIONS = c("missing",
            "Sub-Saharan Africa", 
            "Middle East & North Africa",
            "Latin America & Caribbean",
            "South Asia",
            "East Asia & Pacific",
            "Europe & Central Asia",
            "North America")

