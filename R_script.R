source("R_code/functions.R")

### data preparation
# data input
endo <- read.csv("files/20190726_endophytic_Fungi.csv", h = T,
  sep = "|")

# standardize country names
endo <- get_country(endo)

# assign corrected taxonomy (splitdata set to reduce memory use)
tax <- read.csv("files/taxonomy.csv", h = T, sep = ";")
x <- 0:nrow(tax) %/% 1000
temp <- split(tax, x)
for(i in 1:length(temp)) {
  temp[[i]] <- modify_taxonomy(temp[[i]])
  message(paste("Dataset #", i, " done."))
}
tax <- do.call("rbind", temp)
endo$kingdom <- tax$kingdom
endo$class <- tax$class
endo$order <- tax$order
endo$family <- tax$family
endo$genus <- tax$genus
endo$species <- tax$species
rm(tax, temp)

# remove non-endophyte records
endo <- remove_nonendophytes(endo)

# find out organ of origin of endophytes
endo <- get_plant_part(endo)

# get host data
endo <- correct_host_names(endo)
endo <- get_host_data(endo)

# save data
saveRDS(endo, file = "output/endophytes_data.RData")
endo <- readRDS("output/endophytes_data.RData")

### data analysis
# plot distribution of endophytes records
plot_map(endo, trans = "log10") # output/endophytes_map.png

# plant organs of origin
tab_organ <- table(endo$plant_part)
plot_isolation_organ(tab_organ) # output/isolation_organ.png

# environmental vs. cultures
tab_seq_source <- table(endo$environmental_sample)
plot_seq_source(tab_seq_source) # output/seq_source.png

# taxonomic distribution
plot_tax(endo) # output/endoph_taxa.png, output/endoph_taxa_part.png

# main hosts of endophytes
tab_host <- table(endo$host_family)
plot_host_plant(endo) # output/seq_source.png


### end
