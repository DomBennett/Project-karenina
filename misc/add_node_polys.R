# Dom Bennett
# 20/07/2015
# Add a tip of length 0 to every node in order to make a correlation matrix

library (MoreTreeTools)




tree <- rtree (10)
plot (tree)
# loop through each internal node and add a tip


plot (tree)


nnode <- getSize (tree) + tree$Nnode
intnode <- getSize (tree) + 1
all.node.labels <- paste0 ('n', 1:nnode)
tree$tip.label <- all.node.labels[1:getSize (tree)]
plot (tree)
counter <- 0

# get cormatrix
dmat <- cor (dist.nodes (trees[[i]]))
nnode <- getSize (trees[[i]]) + trees[[i]]$Nnode
all.node.labels <- paste0 ('n', 1:nnode)
colnames (dmat) <- rownames (dmat) <- all.node.labels
tree <- vcv2phylo (lower.tri (dmat))

corobj <- corSymm (1)
corobj <- Initialize.corPhyl
dmat[lower.tri (dmat)], fixed=TRUE
corMatrix (corr=dmat)
summary (corobj)