library(ggplot2)
library(dplyr)

p_DEGs <- read.csv(file = "F:/Bioinformatics/Single Cell Work/Diagram/Heatmap/Final File/common gene logfc.csv")

# SARS Cov-2 DEGs Heatmap
# Custom x and y labels with cexRow and labRow (col respectively)
rnames <- p_DEGs[,1]
mat_data <- data.matrix(p_DEGs[,2:3])
rownames(mat_data) <- rnames

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("Beige", "orange", "gray"))(n = 299)
# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(-1,0,length=100),  # for red
               seq(0.01,0.8,length=100),           # for yellow
               seq(0.81,1,length=100))             # for green
heatmap(mat_data, scale="row", col = my_palette, cexRow = 0.5, Rowv = NA, legend="none")


#mat_data<-as.matrix(All108genes_COVID_ASthema_COPD_CF_IPF_LFCtable)
#mat_data<-as.matrix(p_DEGs)
mat_data
par(mar=c(-1,1))
#pheatmap(mat_data, scale = "row", Rowv = NA, col=my_palette)

#pheatmap
library(pheatmap)
pheatmap(mat_data,  fontsize = 8,treeheight_col = 20)
pheatmap(mat_data,  fontsize = 10,treeheight_col = 20, col=colorRampPalette(c("white","red"))(200))





# Remove "legend" argument from heatmap function
heatmap(mat_data, scale = "row", col = my_palette, cexRow = 0.5, Rowv = NA)

# Correct margin setting
par(mar = c(5, 4, 4, 2) + 0.1)

# Use pheatmap function instead
library(pheatmap)
pheatmap(mat_data, fontsize = 8, treeheight_col = 20)
pheatmap(mat_data, fontsize = 10, treeheight_col = 20, color = colorRampPalette(c("white", "red"))(200))



# Load the pheatmap library
library(pheatmap)

# Define a custom color palette
my_color_palette <- colorRampPalette(c("blue", "white", "red"))(n = 200)

# Create the heatmap with the custom color palette
pheatmap(mat_data, fontsize = 10, treeheight_col = 20, color = my_color_palette)


# Load the pheatmap library
library(pheatmap)

# Define the "viridis" color palette
library(viridis)
my_color_palette <- viridis(200)

# Create the heatmap with the "viridis" color palette
pheatmap(mat_data, fontsize = 10, treeheight_col = 20, color = my_color_palette)



