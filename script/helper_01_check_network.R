library("PlackettLuce")
library("BradleyTerryScalable")
library("igraph")


R <- G[1:length(G),,as.grouped_rankings = FALSE]

R <- as.rankings(R)

adj <- adjacency(R)

adj <- as.vector(adj)

adj <- t(matrix(adj, nrow = ncol(R), ncol = ncol(R)))

dimnames(adj) <- list(dimnames(R)[[2]], dimnames(R)[[2]])

adj <- btdata(adj, return_graph = TRUE)


svg(filename = paste0(output, "network.svg"), width = 12, height = 12, pointsize = 20)
plot.igraph(adj$graph, vertex.size = 5, edge.arrow.size = 0.1) 
dev.off()

