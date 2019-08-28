
# set output folder
if (!dir.exists("output")) dir.create("output", recursive = T)


### retain endophyte records from dataset --------------------------------------
filter_endophytes <- function(data) {
  sel <- c()
  for(i in colnames(data)) {
    sel <- c(sel, grep("endoph", data[, i]))
  }
  sel <- unique(sel)  
  res <- data[sel, ]
  return(res)
}


### lat_lon to separate columns ------------------------------------------------
split_lat_lon <- function(data) {
  coord <- as.character(data$lat_lon)
  lat <- as.numeric(sapply(strsplit(coord, " "), "[", 1))
  y   <- sapply(strsplit(coord, " "), "[", 2)
  lon <- as.numeric(sapply(strsplit(coord, " "), "[", 3))
  x   <- sapply(strsplit(coord, " "), "[", 4)
  y[is.na(y)] <- "-"
  for(i in 1:length(y)) if(y[i] == "S") lat[i] <- -lat[i]
  x[is.na(x)] <- "-"
  for(i in 1:length(x)) if(x[i] == "W") lon[i] <- -lon[i]
  data$lon <- lon
  data$lat <- lat
  return(data)
}


### get standardized country information ---------------------------------------
get_country <- function(data) {
  require(rnaturalearth)

  data <- split_lat_lon(data)
  country <- sub(":.*", "", data$country)
  
  # records with coordinates but no country field
  sel <- which(!is.na(data$lat_lon) & is.na(data$country))
  coord <- na.omit(data[sel, c("lon", "lat")])
  country2 <- coords2country(coord)
  country[sel] <- country2

  # manually change some names
  country[country == "USA"] <- "United States"
  country[country == "UK"] <- "United Kingdom"
  country[country == "Viet Nam"] <- "Vietnam"
  country[country == "Cote d'Ivoire"] <- "Côte d'Ivoire"
  country[country == "Svalbard"] <- "Norway"
  country[country == "Czech Republic"] <- "Czech Rep."
  country[country == "South Korea"] <- "Korea"
  #country[country == "French Guiana"] <- "Guyana"

  # check country names against map names
  map_countries <- ne_countries(scale = "medium", returnclass = "sf")$name
  message("Countries not matching map (removed):")
  print(unique(na.omit(country[!(country %in% map_countries)])))
  
  # add data to file
  data$country2 <- country
  return(data)
}

### converting coordinates to countries ----------------------------------------
# taken from https://stackoverflow.com/questions/41105042/using-coords2country-function-in-r-on-exclusive-economic-zones-not-country-bound
# The single argument to this function, points, is a data.frame in which:
#   - column 1 contains the longitude in degrees
#   - column 2 contains the latitude in degrees
coords2country <- function(points) {
  require(sp)
  require(maps)
  require(rgeos)
  require(maptools)
 
  # prepare a SpatialPolygons object with one poly per country
  countries = map('world', fill = T, col = "transparent", plot = F)
  names = sapply(strsplit(countries$names, ":"), function(x) x[1])

  # clean up polygons that are out of bounds
  filter = countries$x < -180 & !is.na(countries$x)
  countries$x[filter] <- -180

  filter = countries$x > 180 & !is.na(countries$x)
  countries$x[filter] <- 180

  countriesSP <- map2SpatialPolygons(countries, IDs=names,
    #proj4string = CRS("+proj=longlat +datum=wgs84"))
    proj4string = CRS())

  # convert our list of points to a SpatialPoints object
  pointsSP <- SpatialPoints(points,
    proj4string = CRS())
    #proj4string = CRS("+proj=longlat +datum=wgs84"))

  # use 'over' to get indices of the Polygons object containing each point 
  indices <- over(pointsSP, countriesSP)

  # Return the state names of the Polygons object containing each point
  myNames <- sapply(countriesSP@polygons, function(x) x@ID)
  myNames[indices]
}


### plot maps with countries colored by number of records ----------------------
plot_map <- function(data, trans = "sqrt") {
  require(ggplot2)
  require(ggthemes)
  require(rnaturalearth)
  require(sf)

  theme_set(theme_bw())
  #theme_set(theme_map())
  world <- ne_countries(scale = "medium", returnclass = "sf")
 
  # number of records per country
  counts <- table(data$country2)
  records <- rep(0, length(world$name))
  names(records) <- world$name
  for(i in names(counts)[names(counts) %in% names(records)]) {
    records[i] <- counts[i]
  }
  
  # plot map with colored countries
  ggplot(data = world) + geom_sf(aes(fill = records)) +
    scale_fill_distiller(palette = "Blues", direction = 1, trans = trans,
    na.value = "white") + ggtitle("Fungal endophyte records in NCBI GenBank",
    subtitle = paste0(format(sum(!is.na(data$country2)), big.mark = ",",
    scientific = F), " out of ", format(nrow(endo), big.mark = ",",
    scientific = F), " (", round(sum(!is.na(data$country2)) * 100 / nrow(endo),
    1), "%) with country information")) + theme(text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5))
  ggsave(paste0("output/endophytes_map.png"), w = 9, h = 6)
  #ggsave(paste0("output/endophytes_map.pdf"), w = 9, h = 6)
} 



### obtaind data on the plant part of origin of endophytes ---------------------
get_plant_part <- function(data) {
  # terms to differentiate above- and below-grund tissues
  leaf_terms <- c("leaf", "leaves", "frond", "foliar", "peduncle", "twig",
    "foliage", "canopy", "petiole", "needle", "phyllosphere")
  stem_terms <- c("stem", "thallus", "bark", "trunk", "petiole", "branch",
    "phloem", "braches", "cladode", "stalk", "shoot", "moss", "livewort",
    "tiller", "thorn")
  seed_terms <- c("^seed$", "^seeds$", "fruit", "berry", "flower", "chestnuts",
    "ovary", "ovule", "floret")
  root_terms <- c("root", "rhizome", "dark septate endophyte", "DSE", "crown",
    "ectomycorrhiza", "rhizoid", "tuber", "bulb")
  
  # detect organ of origin
  organ <- rep(NA, nrow(data))
  for(i in 1:length(organ)) {    
    leaf <- NA
    for(j in leaf_terms) {
      result <- grep(j, data$isolation_source[i], ignore.case = T)
      if(length(result) != 0) leaf <- "leaf"
      rm(result)
    }
    stem <- NA
    for(j in stem_terms) {
      result <- grep(j, data$isolation_source[i], ignore.case = T)
      if(length(result) != 0) stem <- "stem"
      rm(result)
    }
    seed <- NA
    for(j in seed_terms) {
      result <- grep(j, data$isolation_source[i], ignore.case = T)
      if(length(result) != 0) seed <- "seed, fruit, or flower"
      rm(result)
    }
    root <- NA
    for(j in root_terms) {
      result <- grep(j, data$isolation_source[i], ignore.case = T)
      if(length(result) != 0) root <- "root"
      rm(result)
    }
    organ[i] <- paste(leaf, stem, seed, root, sep = "; ")
  }
  organ <- gsub("; NA", "", organ)
  organ <- gsub("NA; ", "", organ)
  organ[grep("; root", organ)] <- "systemic"
  organ[grep("; ", organ)] <- "above-ground organs"
  data$plant_part <- organ

  return(data)
}


### modify names in taxonomy file ----------------------------------------------
modify_taxonomy <- function(x) {
  ncat <- ncol(x)
  for(i in rownames(x)) {
    for(j in 2:(ncat - 1)) {
      if(any(grep("unclassified", x[i, j])) &
        any(grep("unclassified", x[i, j - 1], invert = T))) {
        tag <- paste0("unclassified_", x[i, j - 1])
        if(!(tag %in% levels(x[, j]))) {
          levels(x[, j]) <- c(levels(x[, j]), tag)
          x[i, j] <- tag
        } else {
          x[i, j] <- tag
        }
      } else if(any(grep("unclassified", x[i, j])) &
        any(grep("unclassified", x[i, j - 1]))) {
        tag <- as.character(x[i, j - 1])
        if(!(tag %in% levels(x[, j]))) {
          levels(x[, j]) <- c(levels(x[, j]), tag)
          x[i, j] <- tag
        } else {
          x[i, j] <- tag
        }      
      }
    }
    if(any(grep("unclassified", x[i, ncat - 1]))) {
      tag <- as.character(x[i, ncat - 1])
      if(!(tag %in% levels(x[, ncat]))) {
        levels(x[, ncat]) <- c(levels(x[, ncat]), tag)
        x[i, ncat] <- tag
      } else {
        x[i, ncat] <- tag
      }
    } else if(any(grep("_sp", x[i, ncat], invert = T))) {
      tag <- paste0(x[i, ncat - 1], "_sp")
      if(!(tag %in% levels(x[, ncat]))) {
        levels(x[, ncat]) <- c(levels(x[, ncat]), tag)
        x[i, ncat] <- tag
      } 
    } 
    if(any(grep("unclassified", x[i, ncat])) &
      any(grep("unclassified", x[i, ncat - 1], invert = T))) {
        tag <- paste0(x[i, ncat - 1], "_sp")
        if(!(tag %in% levels(x[, ncat]))) {
          levels(x[, ncat]) <- c(levels(x[, ncat]), tag)
          x[i, ncat] <- tag
        } else {
          x[i, ncat] <- tag
        }
    }
  }
  return(x)
}


### retrieve host data ---------------------------------------------------------
get_host_data <- function(data) {
  require(Taxonstand)
  host_list <- unique(data$host[!is.na(data$host)])
  host_list_corrected <- TPL(host_list)
  data$host_species <- data$host
  data$host_family <- rep(NA, nrow(data))
  for(i in 1:length(host_list)) {
    data$host_species[data$host_species == host_list[i]] <-
      paste(host_list_corrected$New.Genus[i],
      host_list_corrected$New.Species[i])
    data$host_family[data$host_species == host_list[i]] <-
      host_list_corrected$Family[i]
  }
  data$host_species <- gsub(" NA$", " sp.", data$host_species)
  return(data)
}

### manually correct host names ------------------------------------------------
correct_host_names <- function(data) {
  data$host_original <- data$host

  # remove 'sp.' epithet
  data$host <- gsub(" sp.", "", data$host, fixed = T)
  data$host <- trimws(data$host)
  data <- droplevels(data)
  
  # manually correct names
  names <- data.frame(
    old = c(
      "soybean",
      "Norway spruce",
      "tobacco",
      "soybean (Conquista)",
      "Scots pine",
      "soybean (M-Soy 6101)",
      "isolated from twigs of Quercus brantii",
      "Pinus tabulaeformis",
      "Persian oak",
      "sugarcane cultivar IMI-1",
      "mango",
      "maize",
      "patchouli",
      "silver birch",
      "rice",
      "wild rice",
      "sugarcane cultivar SP80-1842",
      "palm trees",
      "Norway spruce seedling",
      "pearl millet",
      "tomato",
      "olive cultivar Cobrancosa",
      "soybean; Conquista (MG/BR46)",
      "coffee; Catigua MG2",
      "coffee; Catuai Vermelho IAC 144",
      "soybean; Vencedora (BRSMG 68)",
      "peanut",
      "lotus leaves",
      "poplar",
      "banana",
      "wild Oryza sativa",
      "sugar cane",
      "sugar beet",
      "Salix (willow)",
      "Citrus (pomelo)",
      "mung bean",
      "apple tree",
      "Taraxacum (herb)",
      "rubber tree",
      "tobacco plant",
      "finger millet",
      "coffee arabica",
      "NMiscanthus floridulus",
      "Japanese cedar",
      "milk thistle",
      "acid lime",
      "Acid lime",
      "molokia plant",
      "Chinese maple",
      "Quercus phillyraeoides",
      "white oak"
    ), new = c(
      "Glycine max",
      "Picea abies",
      "Nicotiana",
      "Glycine max",
      "Pinus sylvestris",
      "Glycine max",
      "Quercus brantii",
      "Pinus tabuliformis",
      "Quercus brantii",
      "Saccharum officinarum",
      "Mangifera indica",
      "Zea mays",
      "Pogostemon cablin",
      "Betula pendula",
      "Oryza sativa",
      "Zizania",
      "Saccharum officinarum",
      "Phoenix dactylifera",
      "Picea abies",
      "Pennisetum glaucum",
      "Solanum lycopersicum",
      "Olea europaea",
      "Glycine max",
      "Coffea arabica",
      "Coffea arabica",
      "Glycine max",
      "Arachis hypogaea",
      "Lotus",
      "Populus",
      "Musa",
      "Oryza sativa",
      "Saccharum officinarum",
      "Beta vulgaris",
      "Salix",
      "Citrus",
      "Vigna radiata",
      "Malus domestica",
      "Taraxacum",
      "Hevea brasiliensis",
      "Nicotiana",
      "Eleusine coracana",
      "Coffea arabica",
      "Miscanthus floridulus",
      "Cryptomeria japonica",
      "Silybum marianum",
      "Citrus hystrix",
      "Citrus hystrix",
      "Corchorus olitorius",
      "Acer palmatum",
      "Quercus phillyreoides",
      "Quercus alba")
  )
  for(i in 1:nrow(names)) {
    #if(!(names[i, "new"] %in% levels(data$host))) {
    #  levels(data$host) <- c(levels(data$host), names[i, "new"])
    #}
    data$host[data$host == as.character(names[i, "old"])] <- names[i, "new"]
  }
  data <- droplevels(data)

  return(data)
}

### remove non endophyte records -----------------------------------------------
remove_nonendophytes <- function(data) {
  terms <- c("soil", "dung", "hindgut", "gut of larva", "gallery of beetle",
    "gallery", "exoskeleton of a moth (Lacinipolia sp.) visiting Paspalum dilatatum inflorescences with evidence of infection by Claviceps paspali",
    "cornea scrape from woman with right eye scratched by wood",
    "corneal ulcer", "contaminate of fungal culture", "coral", "compost",
    "clinical sample (eye)", "abscess", "surrounding air sample", "sponge",
    "soil in maize field","pulp mill slime", "mushroom", "mushroom casing",
    "indoor air", "human cornea", "Homo sapiens", "Galleria mellonella",
    "oral secretions in adult Dendroctonus rufipennis",
    "surface-sterilized lichen thallus",
    "recovered in culture on 2% malt extract agar from surface-sterilized lichen thallus",
    "endolichenic fungus recovered in culture on 2% malt extract agar from surface-sterilized lichen thallus",
    "lichen growing on Acer pseudoplatanus", "sterile lichen")
  for(i in terms) {
    data <- data[!(data$isolation_source %in% i), ]
    data <- data[!(data$host %in% i), ]
  }
  data <- droplevels(data)
  return(data)
}


### function to generate distinct colors ---------------------------------------
# adapted from https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
color <- function(n, start = 1, random = F, last_gray = T) {
  require(RColorBrewer)
  
  # 433 colors
  # x <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  
  # 74 colors
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  x <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
    rownames(qual_col_pals)))

  result <- x[start:(start + n - 1)]
  if(random) result <- sample(x, n)
  if(last_gray) result[length(result)] <- gray(0.75)
  
  return(result)  
}


### plot pie chart with organs of origin of endophytes -------------------------
plot_isolation_organ <- function(data) {
  names(data)[which(names(data) == "NA")] <- "unknown"
  data <- sort(data, decreasing = T)
  data <- data[c(3, 4, 5, 6, 7, 2, 1)]
  pie_col <- c("chartreuse4", "chartreuse3", "chartreuse2",
    "chartreuse1", "brown1", "darkgoldenrod3", gray(0.7))
  pie_lab <- paste0(round(data * 100 / sum(data), 1), " %")
  pie_lab[c(4, 5)] <- rep("", 2)
  png("output/isolation_organ.png", h = 3, w = 5, units = "in", res = 600,
    pointsize = 12)
  par(mar = rep(2, 4), mfrow = c(1, 2), xpd = T)
  pie(data, labels = pie_lab, col = pie_col, border = "white",
    clockwise = T)
  mtext("Plant organ of origin", side = 3, line = -1)
  par(mar = rep(1, 4))
  plot(1, 1, type = "n", axes = F)
  legend("center", legend = names(data), bty = "n", fill = pie_col,
    border = NA, y.intersp = 0.75)
  dev.off()  
}


### plot pie chart with sequence source ----------------------------------------
plot_seq_source <- function(data) {
  names(data) <- c("specimen", "environmental")
  pie_col <- color(2, last_gray = F)
  pie_lab <- paste0(round(data * 100 / sum(data), 1), " %")
  png("output/seq_source.png", h = 3, w = 5, units = "in", res = 600,
    pointsize = 12)
  par(mar = rep(2, 4), mfrow = c(1, 2), xpd = T)
  pie(data, labels = pie_lab, col = pie_col, border = "white",
    clockwise = T)
  mtext("Source of sequences", side = 3, line = -1)
  par(mar = rep(0, 4))
  plot(1, 1, type = "n", axes = F)
  legend("left", legend = names(data), bty = "n", fill = pie_col,
    border = NA, y.intersp = 0.75)
  dev.off()  
}


### plot bar chart with taxonomy -----------------------------------------------
plot_tax <- function(data, n = 15) {
  tab_tax <- table(data$order)
  tab_tax_part <- table(data$order, data$plant_part)
  tab_tax <- sort(tab_tax, decreasing = T)
  others <- sum(tail(tab_tax, length(tab_tax) - n))
  names(others) <- "others"
  others_part <- tab_tax_part[names(tail(tab_tax, length(tab_tax) - n)), ]
  tab_tax <- head(tab_tax, n)
  tab_tax_part <- tab_tax_part[names(tab_tax), ]
  tab_tax <- c(tab_tax, others)
  tab_tax_part <- rbind(tab_tax_part, others = colSums(others_part))

  # plot overall frequency
  png("output/endoph_taxa.png", h = 4, w = 8, units = "in", res = 600,
    pointsize = 12)
  par(mar = c(8, 6, 2, 2), font.main = 1)
  x <- barplot(tab_tax, names.arg = NA, border = F, las = 2, col = gray(0.5),
    main = "Main endophyte orders in NCBI GenBank")
  text(x, y = -1000, gsub("_", " ", names(tab_tax)), xpd = TRUE, srt = 45,
    adj = 1)
  mtext("records (n)", side = 2, line = 4)
  dev.off()

  # plot frequency per plant part
  colnames(tab_tax_part)[which(colnames(tab_tax_part) == "NA")] <- "unknown"
  tab_tax_above <- rowSums(tab_tax_part[, c(1, 2, 5, 6)])
  tab_tax_root <- tab_tax_part[, "root"]
  tab_tax_others <- rowSums(tab_tax_part[, c("systemic", "unknown")])
  tab_tax_part <- rbind(tab_tax_above, tab_tax_root, tab_tax_others)
  rownames(tab_tax_part) <- c("above-ground", "below-ground",
    "systemic/unknown")
  png("output/endoph_taxa_part.png", h = 4, w = 8, units = "in", res = 600,
    pointsize = 12)
  par(mar = c(8, 6, 2, 2), font.main = 1)
  part_col <- c("chartreuse4", "darkgoldenrod3", gray(0.7))
  par(mar = c(8, 6, 2, 2), font.main = 1)
  x <- barplot(tab_tax_part, names.arg = rep("", ncol(tab_tax_part)),
    border = F, las = 2, col = part_col,
    main = "Main endophyte orders in NCBI GenBank")
  text(x, y = -1000, gsub("_", " ", names(tab_tax)), xpd = TRUE, srt = 45,
    adj = 1)
  mtext("records (n)", side = 2, line = 4)
  legend("top", legend = rownames(tab_tax_part), fill = part_col,
    border = F, bty = "n", y.intersp = 0.75)  
  dev.off()
}


### plot bar chart with hosts --------------------------------------------------
plot_host_plant <- function(data, n = 15) {
  require(plotrix)
  data <- table(data$host_family)
  data <- sort(data, decreasing = T)
  unknown <- data[which(names(data) == "")]
  names(unknown) <- "unknown"
  data <- data[-which(names(data) == "")]
  others <- sum(tail(data, length(data) - n))
  names(others) <- "others"
  data <- head(data, n)
  data <- c(data, others, unknown)

  # plot overall frequency
  png("output/endoph_hosts.png", h = 4, w = 8, units = "in", res = 600,
    pointsize = 12)
  par(mar = c(8, 6, 2, 2), font.main = 1)
  x <- barplot(data, names.arg = NA, border = F, las = 2, col = gray(0.5),
    main = "Main endophyte hosts in NCBI GenBank")
  text(x, y = -1000, gsub("_", " ", names(data)), xpd = TRUE, srt = 45,
    adj = 1)
  mtext("records (n)", side = 2, line = 4)
  dev.off()
}


### end




