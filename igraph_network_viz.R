#Packages
library(igraph)
library(tidyr)
library(ggplot2)
getwd()

#NETWORK STATIC VISUALIZATIONS-------------------------------------------------------------------------------------

#Extract names of different networks from the folder containing the nodes and
#edges (as dataframes) of the communities extracted with the 'KinaseNetworks'
#jupyter-notebook
network_names <- grep('edges', list.files(path = 'output_dataframes'), value=T)
network_names <- gsub('_edges\\.csv', '', network_names)

#Loop through each network
for (network in network_names[3]){

  #Load edges data
  edges <- read.csv( paste0('output_dataframes/', network, '_edges.csv'), row.names = 1 )
  #Load nodes data
  nodes <- read.csv( paste0('output_dataframes/', network, '_nodes.csv'), row.names = 1 )
  #Build network
  net <- graph.data.frame(edges, nodes, directed=F)
  #Network settings
  E(net)$width <- E(net)$weight/15
  V(net)$label.cex = 0.5
  V(net)$size <- 14
  if(network == "AKT.AZD5363_vs_Control"| network=="PI3K.GDC0941_vs_Control"){
  nodes_highlight <- match( c('PIK3CA', 'AKT1_2', 'MTOR', 'RPS6KB1'), nodes$id )
  nodes_highlight <- nodes_highlight[!is.na(nodes_highlight)]
  l_factor <- 1.5
  col <- "#C5E5E7"
  } else if (network=="ERK.GDC0994_vs_Control"| network=="MEK.Trametinib_vs_Control" ){
  nodes_highlight <- match( c('MAP2K1', 'MAPK1_3'), nodes$id )
  nodes_highlight <- nodes_highlight[!is.na(nodes_highlight)]
  col <- c("#C5E5E7")
  l_factor <- 0.6
  } else {
  AKT <- match( c('PIK3CA', 'AKT1_2', 'MTOR', 'RPS6KB1'), nodes$id )
  AKT <- AKT[!is.na(AKT)]
  MAP2K1 <- match( c('MAP2K1', 'TNK2'), nodes$id )
  MAP2K1 <- MAP2K1[!is.na(MAP2K1)]
  nodes_highlight <- list(AKT, MAP2K1)
  col <- c("#C5E5E7","#ECD89A")
  E(net)$width <- E(net)$weight*3
  l_factor <- 1
  }
  l <- layout.fruchterman.reingold(net)
  
  #Export network
  pdf(file= paste0("output_visualizations/", network ,"_networkviz_static.pdf"))
  plot(net, vertex.label.family='Helvetica', edge.arrow.size=.4, rescale=F, layout=l*l_factor, mark.groups =nodes_highlight, mark.col = col, mark.border = NA )
  dev.off()
}

#Interaction plot - save manually as it 
plot(net, vertex.label.family='Helvetica', edge.arrow.size=.4, layout=l, mark.groups =nodes_highlight, mark.col = col, mark.border = NA )

#GROUPED HORIZONTAL BARPLOTS---------------------------------------------------------------------------------------------------------------------------------------

#Read edges dataset
edges_df <- read.csv('PROJECT_DATASET_2.csv')
#Convert positive edges to NA
edges_df[,2:ncol(edges_df)][edges_df[,2:ncol(edges_df)]>0] <- NA
#Convert negative values to positive
edges_df[,2:ncol(edges_df)] <- abs(edges_df[,2:ncol(edges_df)])
#Separate first column by dot
edges_df <- edges_df %>% separate(X, sep = "\\.", into = c('from','to'), remove = T)

#PIK3CA AKT
sel <- c('AKT1_2','MTOR','PIK3CA','RPS6KB1')
edges_akt <- subset(edges_df, from %in% sel & to %in% sel)[,c(1,2,4,5)]
edges_akt_long <- pivot_longer(edges_akt, colnames(edges_akt)[c(3,4)])
edges_akt_long$`absolute z-score` <- edges_akt_long$value
edges_akt_long$edge <- paste0(edges_akt_long$from, '-', edges_akt_long$to)

# Change the colors manually
p <- ggplot(data=edges_akt_long, aes(x=edge, y=`absolute z-score`, fill=name)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()
# Use brewer color palettes
p + scale_fill_brewer(palette="Blues")

plot(1:10, main = expression("GDC0941"["(PIK3CA)"]^"-"))

#MAP2K1 MAPK1_3



#--------------------------------------------------------------------------------------------------------------------------


