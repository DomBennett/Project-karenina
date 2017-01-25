

library(MoreTreeTools)

plot(tree, show.tip.label = FALSE)

tree <- rtree(1000)
eds <- calcED(tree)
tip <- rownames(eds)[eds == min(eds[,1])]
edges <- getEdges(tree, tips=tip, type=2)
for(edge in edges) {
  es <- getEdges(tree, node=tree$edge[edge, 2])
  tree$edge.length[es] <- tree$edge.length[es]*2
}
plot(tree, show.tip.label = FALSE, no.margin=TRUE)
