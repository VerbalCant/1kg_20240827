# This script is a modified version of the script found in the excellent
# Biostars tutorial by Kevin Blighe on building a PCA bi-plot at https://www.biostars.org/p/335605/.
# The only requirement I've added is dplyr because life without dplyr is not worth living.
# Well, I also added the dependency on the IGSR populations table, but it's in the repo.
# 2024-08-26 Alaina Hardie

options(scipen=100, digits=3)
library(dplyr)

# To produce these tables referenced below, follow the workflow instructions in the tutorial linked above.
# I have also committed my own results to this repo, so this R script should run fine with the local versions
# if you have cloned the whole repo.
eigenvec <- read.table('plink.eigenvec', header = FALSE, skip=0, sep = ' ')
rownames(eigenvec) <- eigenvec[,2]
eigenvec <- eigenvec[,3:ncol(eigenvec)]
colnames(eigenvec) <- paste('Principal Component ', c(1:20), sep = '')
special_populations <- c("PEL", "CLM", "MXL", "PUR")

PED <- read.table('20130606_g1k.ped', header = TRUE, skip = 0, sep = '\t')
PED <- PED[which(PED$Individual.ID %in% rownames(eigenvec)), ]
PED <- PED[match(rownames(eigenvec), PED$Individual.ID),]
PED <- PED %>%
  group_by(Population) %>%
  group_modify(~ {
    if (.y$Population %in% special_populations) {
      .x  # Keep all individuals for these populations
    } else {
      slice_sample(.x, n = min(100, nrow(.x)))  # Limit to 100 for others
    }
  }) %>%
  ungroup()

eigenvec <- eigenvec[rownames(eigenvec) %in% PED$Individual.ID, ]

eigenvec <- eigenvec[match(PED$Individual.ID, rownames(eigenvec)), ]

if (!all(PED$Individual.ID == rownames(eigenvec))) {
  stop("PED Individual IDs do not match eigenvec rownames.")
}

require('RColorBrewer')

# Load this as a reference so I don't have to type out colours and titles.
igsr_data <- read.csv("igsr_populations.tsv", sep="\t", stringsAsFactors = FALSE)

# Choose all of the populations with population codes
popnames = c("FIN", "CHS", "KHV", "BEB", "PUR", "ACB", "GWW", "ASW", "YRI", "GWD", "JPT", "MSL", "CEU", "ESN", "CHB", "CLM", "CDX", "PEL", "PJL", "IBS", "TSI", "MXL", "LWK", "GIH", "GWF", "STU", "ITU", "GWJ", "GBR", "MKK")
filtered_data <- igsr_data[igsr_data$Population.code %in% popnames, ]
poptitles <- paste(filtered_data$Population.code, filtered_data$Population.description, sep = ": ")

custom_colors <- c("PUR" = "red",  
                   "CLM" = "green", 
                   "PEL" = "blue",
                   "MXL" = "orange")

colournames <- filtered_data$Superpopulation.display.colour
names(colournames) <- filtered_data$Population.code

for (pop in names(custom_colors)) {
  colournames[pop] <- custom_colors[pop]
}

colournames <- unname(colournames)

# from: http://www.internationalgenome.org/category/population/
PED <- PED[PED$Population %in% popnames, ]
col <- colournames[as.numeric(factor(PED$Population, levels = popnames))]

project.pca <- eigenvec
summary(project.pca)

par(mar = c(5,5,5,5), cex = 2.0,
    cex.main = 7, cex.axis = 2.75, cex.lab = 2.75, mfrow = c(1,2))

# Write it to a PDF locally.
pdf("PCA_Plot_1_vs_2.pdf", width = 16, height = 12)

# Define custom shapes for our Possibly Proximal Populations of interest
custom_shapes <- c("PUR" = 17,  
                   "CLM" = 18,  
                   "PEL" = 23,  
                   "MXL" = 15)  

# Default shape for all other populations
default_shape <- 20

shapes <- rep(default_shape, nrow(PED))

# Assign custom shapes to specific populations
for (pop in names(custom_shapes)) {
  shapes[PED$Population == pop] <- custom_shapes[pop]
}
point_sizes <- rep(1, nrow(PED)) 
point_sizes[PED$Population == "PEL"] <- 2

pdf("1kg_PCA_Plot_PC1_vs_PC2.pdf", width = 16, height = 12)

plot(project.pca[,1], project.pca[,2],
     type = 'n',
     main = 'Selected 1kg populations',
     adj = 0.5,
     xlab = 'PC1',
     ylab = 'PC2',
     font = 2,
     font.lab = 2)

points(project.pca[,1], project.pca[,2], col = col, pch = shapes, cex = point_sizes)

legend_data <- data.frame(poptitles = poptitles, colournames = colournames)
legend_data <- legend_data[order(legend_data$colournames), ]

legend('bottom', 
       bty = 'n', 
       cex = 1, 
       ncol = 2, 
       title = '', 
       legend = legend_data$poptitles, 
       fill = legend_data$colournames, 
       xpd = TRUE)

par(mar = c(5, 5, 4, 10))
dev.off() # close after PDF save