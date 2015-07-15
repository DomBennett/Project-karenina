## No copyright, no warranty
## Dominic John Bennett
## Functions for calculating LFI
## 07/02/2014

## Dependencies
require (ape)
require (geiger)
require (plyr)
library (graphics)

## Functions
getTimeslice <- function (phylo, time.slice) {
  # identify all nodes and edges that need to be kept
  # all edges that have 1 node on older than time.slice
  keep.nodes <- which (phylo$node.age >= time.slice)
  keep.edges <- which (phylo$edge[ ,1] %in% keep.nodes)
  # remove any tip nodes that don't pass through timeslice
  tip.edges <- which (phylo$edge[keep.edges, 2] <= length (phylo$tip.label))
  tip.node.ages <- phylo$node.age[phylo$edge[keep.edges[tip.edges], 2]]
  drop.tip.edges <- tip.edges[tip.node.ages > time.slice]
  if (length (drop.tip.edges) > 0) {
    keep.edges <- keep.edges[-drop.tip.edges]
  }
  # create an edge matrix
  edge.matrix <- phylo$edge[keep.edges, ]
  edge.lengths <- phylo$edge.length[keep.edges]
  edge.changes <- phylo$edge.changes[keep.edges]
  tips <- edge.matrix[!edge.matrix[ ,2] %in% edge.matrix[ ,1],2]
  tip.labels <- phylo$node.label[tips]
  n.node <- as.integer (length (unique (edge.matrix[ ,1])))
  # relabel all nodes
  new.edge.matrix <- edge.matrix
  new.edge.matrix[edge.matrix[ ,2] %in% tips,2] <- 1:length (tips)
  new.edge.matrix[!new.edge.matrix[ ,1] %in% new.edge.matrix[ ,2],1] <- length (tips) + 1
  old.nodes <- unique (new.edge.matrix[new.edge.matrix[ ,1] > length (tips) + 1,1])
  new.nodes <- (length (tips) + 2):(length (old.nodes) + (length (tips) + 1))
  for (i in 1:length (old.nodes)) {
    new.edge.matrix[new.edge.matrix[ ,1] == old.nodes[i],1] <- new.nodes[i]
    new.edge.matrix[new.edge.matrix[ ,2] == old.nodes[i],2] <- new.nodes[i]
  }
  storage.mode(new.edge.matrix) <- "integer"
  # create new phylo object
  new.phylo <- list (edge = new.edge.matrix, tip.label = tip.labels,
                     edge.length = edge.lengths, Nnode = n.node)
  class (new.phylo) <- 'phylo'
  new.phylo$edge.changes <- edge.changes
  # identify corresponding nodes
  old.nodes <- unique (c (edge.matrix[ ,1], edge.matrix[ ,2]))
  new.nodes <- unique (c (new.edge.matrix[ ,1], new.edge.matrix[ ,2]))
  node.labels <- phylo$node.label[old.nodes]
  node.ages <- phylo$node.age[old.nodes]
  new.phylo$node.label <- node.labels[order(new.nodes)]
  new.phylo$node.age <- node.ages[order(new.nodes)]
  # identify all tip nodes and re-lengthen
  for (i in 1:length (new.phylo$tip.label)) {
    diff <- time.slice - new.phylo$node.age[i]
    edge <- which (new.phylo$edge[ ,2] == i)
    edge.length <- new.phylo$edge.length[edge]
    new.edge.length <- edge.length - diff
    new.phylo$edge.length[edge] <- new.edge.length
    edge.change <- new.phylo$edge.changes[edge]
    new.edge.change <- edge.change * new.edge.length/edge.length
    new.phylo$edge.changes[edge] <- new.edge.change
    new.phylo$node.age[i] <- time.slice
  }
  new.phylo
}

safeFromJSON <- function (url, max.trys = 10) {
  # Wrapper for fromJSON
  trys <- 0
  while (trys < max.trys) {
    json.obj <- try (fromJSON(url)[[1]], silent = TRUE)
    if (class(json.obj) == 'try-error') {
      cat ('---- Connection failed: trying again ----\n')
      trys <- trys + 1
      Sys.sleep (10)
    } else {
      return (json.obj)
    }
  }
  stop ("Failed to connect, server may be down.")
}

findParent <- function (name) {
  # Pull lineage for a fossil record
  extractRankNumber <- function (rec) {
    if (rec$sta == "belongs to") {
      return (rec$rnk)
    }
  }
  extractNames <- function (rec) {
    if (rec$sta == "belongs to") {
      return (rec$nam)
    }
  }
  url <- "http://paleobiodb.org/data1.1/taxa/"
  url <- paste0(url, "list.json?name=", name, "&rel=all_parents")
  json.obj <- safeFromJSON (url)
  if (length (json.obj) > 0) {
    return (rev (ldply (.data = json.obj, .fun = extractNames)[ ,1])) # returned in order
  } else {
    return (c ())
  }
}

findClade <- function (lineages) {
  for (i in length (lineages[[1]]):1) {
    subj <- lineages[[1]][i]
    j <- 2
    while (TRUE) {
      if (j > length (lineages)) {
        success <- TRUE
        break
      }
      query <- lineages[[j]]
      if (subj %in% query) {
        j <- j + 1
      } else {
        success <- FALSE
        break
      }
    }
    if (success) {
      return (subj)
    }
  }
  NA
}

labelNodes <- function (phylo) {
  # Use GNR to label all nodes in a phylogeny
  taxa.res <- taxaResolve (phylo$tip.label)
  nodes <- 1:(length (phylo$tip.label) + phylo$Nnode)
  node.label <- rep (FALSE, length (nodes))
  node.label[1:length (phylo$tip.label)] <-
    unlist(lapply (strsplit(phylo$tip.label, " "), function (x) x[1]))
  for (i in (length (phylo$tip.label) + 1):length (nodes)) {
    descendants <- nodeDescendants(phylo, node = i)
    genus.names <- unlist(lapply (strsplit(descendants, " "), function (x) x[1]))
    if (all (genus.names == genus.names[1])) {
      node.label[i] <- genus.names[1]
    } else {
      lineages <- as.character (taxa.res[taxa.res$search.name %in% descendants, "lineage"])
      lineages <- strsplit (lineages, "\\|")
      lineages <- lineages[!is.na (lineages)]
      if (length (lineages) > 0) {
        node.label[i] <- findClade (lineages)
      }
    }
  }
  phylo$node.label <- node.label
  phylo
}

palaeoPull <- function (clade, limit = 'all') {
  # Pull all records for a clade name PBDB
  url <- paste0 ("http://paleobiodb.org/data1.1/occs/list.json?base_name=",
                 clade,"&limit=", limit)
  json.obj <- safeFromJSON(url)
  extractAgeName <- function (rec) {
    data.frame (name = rec$tna, max.age = rec$eag, min.age = rec$lag)
  }
  downloaded <- ldply (.data = json.obj, .fun = extractAgeName)
  downloaded <- downloaded[!duplicated (downloaded), ]
  downloaded
}

addFossilsToPhylogeny <- function (phylo, records, ex.age = 'min.age') {
  # Add PBDB fossil records to a phylogeny by name and age matching
  #
  # Args
  #  phylo: phylogeny on which fossils are to be added
  #  records: fossil records from PBDB (from palaeoPull)
  #  ex.age: estimate of extinction, either min.age or max.age
  #
  # Return
  #  phylo
  phylo.env <- new.env ()
  local (findEdge <- function (name, age) {
    # Finding matching node in phylogeny by name, then by age
    # Add to oldest matched named node
    matching.nodes <- which (phylo$node.label == name)
    matching.node.ages <- phylo$node.age[matching.nodes]
    node.1 <- matching.nodes[matching.node.ages == max (matching.node.ages)]
    if (length (node.1) > 1) {
      # if more than one matching named node of same age, choose at random
      node.1 <- sample (node.1, 1)
    }
    if (sum (matching.node.ages > age) > 0) {
      node.2 <- phylo$edge[phylo$edge[ ,1] == node.1, 2]
      if (sum (node.2 %in% matching.nodes) != 1) {
        # if all or no decendants of node.1 match name, choose at random
        node.2 <- sample (node.2, 1)
      } else {
        node.2 <- node.2[node.2 %in% matching.nodes]
      }
    } else {
      while (TRUE) {
        node.age <- phylo$node.age[node.1]
        if (node.age > age) {
          break
        }
        node.2 <- node.1
        node.1 <- phylo$edge[phylo$edge[ ,2] == node.1, 1]
      }
    }
    which (phylo$edge[, 1] == node.1 & phylo$edge[ ,2] == node.2)
  }, env = phylo.env)
  local (matchFossilToPhylogeny <- function (record) {
    name <- as.character(record["name"])
    age <- as.numeric (record[ex.age])
    if (name %in% phylo$tip.label) {
      return (NA)
    }
    genus.name <- strsplit(as.character (name), " ")[[1]][1]
    if (genus.name %in% phylo$node.label) {
      lineage <- genus.name
      edge <- findEdge (genus.name, age)
      return (data.frame (name, age, lineage, edge))
    }
    parent <- findParent(genus.name)
    if (length (parent) == 0) {
      print (paste0 ("No lineage record for [", name, "] ..."))
      return (NA)
    }
    matching.clade <- parent[match (TRUE, parent %in% phylo$node.label)]
    if (is.na (matching.clade)) {
      return (NA)
    }
    lineage <- paste (parent[1:which (parent == matching.clade)], collapse = "|")
    edge <- findEdge (matching.clade, age)
    data.frame (name, age, lineage, edge)
  }, env = phylo.env)
  addFossil <- local (function (i) {
    record <- records[i, ]
    record <- matchFossilToPhylogeny (record)
    #print (record)
    if (!is.na (record[[1]])) {
      phylo.edge <- as.integer (record['edge'])
      tip.name <- as.character (record[['name']])
      tip.age <- as.numeric (record['age'])
      if (tip.age == 0) {
        tip.age <- 0.01
      }
      node.age <- phylo$node.age[phylo$edge[phylo.edge, 1]]
      node.age <- node.age - 0.01
      node.label <- as.character (record[['lineage']])
      phylo <<- addTip (phylo, phylo.edge, tip.name, tip.age, node.age, node.label)
    }
  }, env = phylo.env)
  dups <- records$name[duplicated (records$name)]
  for (dup in dups) {
    max.age <- max (records[records$name %in% dup, 'max.age'])
    min.age <- min (records[records$name %in% dup, 'min.age'])
    record <- data.frame (name = dup, min.age, max.age)
    records <- rbind (records[!records$name %in% dup, ], record)
  }
  phylo.age <- max (phylo$node.age)
  records <- records[records[ex.age] < phylo.age, ]
  records <- records [order (records[ex.age], decreasing = TRUE), ]
  local (records <- records, env = phylo.env)
  local (phylo <- phylo, env = phylo.env)
  m_ply (.data = data.frame (i = 1:nrow (records)),
         .fun = addFossil, .expand = FALSE,
         .progress = create_progress_bar (name = "time"))
  get (x = 'phylo', envir = phylo.env)
}
