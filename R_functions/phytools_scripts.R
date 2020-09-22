force.ultrametric=function (tree, method = c("nnls", "extend")) 
{
    method <- method[1]
    if (method == "nnls") 
        tree <- nnls.tree(cophenetic(tree), tree, rooted = TRUE, 
            trace = 0)
    else if (method == "extend") {
        h <- diag(vcv(tree))
        d <- max(h) - h
        ii <- sapply(1:Ntip(tree), function(x, y) which(y == 
            x), y = tree$edge[, 2])
        tree$edge.length[ii] <- tree$edge.length[ii] + d
    }
    else cat("method not recognized: returning input tree\n\n")
    tree
}

nnls.tree = function (dm, tree, rooted = FALSE, trace = 1, weight = NULL, 
    balanced = FALSE) 
{
    if (is.rooted(tree) & rooted == FALSE) {
        tree <- unroot(tree)
        warning("tree was rooted, I unrooted the tree!")
    }
    tree <- reorder(tree, "postorder")
    if (balanced) {
        if (!is.binary(tree)) 
            stop("tree must be binary")
        weight <- rowSums(designTree(unroot(tree)))
    }
    dm <- as.matrix(dm)
    k <- dim(dm)[1]
    labels <- tree$tip.label
    dm <- dm[labels, labels]
    y <- dm[lower.tri(dm)]
    if (rooted) 
        X <- designUltra(tree)
    else X <- designUnrooted2(tree)
    if (!is.null(weight)) {
        y <- y * sqrt(weight)
        X <- X * sqrt(weight)
    }
    lab <- attr(X, "nodes")
    if (any(is.na(y))) {
        ind <- which(is.na(y))
        X <- X[-ind, , drop = FALSE]
        y <- y[-ind]
    }
    Dmat <- crossprod(X)
    dvec <- crossprod(X, y)
    betahat <- as.vector(solve(Dmat, dvec))
    betahattmp <- betahat
    bhat <- numeric(max(tree$edge))
    bhat[as.integer(lab)] <- betahat
    betahat <- bhat[tree$edge[, 1]] - bhat[tree$edge[, 2]]
    if (!any(betahat < 0)) {
        RSS <- sum((y - (X %*% betahattmp))^2)
        if (trace) 
            print(paste("RSS:", RSS))
        attr(tree, "RSS") <- RSS
        tree$edge.length <- betahat
        return(tree)
    }
    n <- dim(X)[2]
    l <- nrow(tree$edge)
    lab <- attr(X, "nodes")
    ind1 <- match(tree$edge[, 1], lab)
    ind2 <- match(tree$edge[, 2], lab)
    Amat <- matrix(0, 2, l)
    Amat[1, ] <- 1
    Amat[2, ] <- -1
    Aind <- matrix(0L, 3, l)
    Aind[1, ] <- 2L
    Aind[2, ] <- as.integer(ind1)
    Aind[3, ] <- as.integer(ind2)
    if (any(is.na(Aind))) {
        na_ind <- which(is.na(Aind), arr.ind = TRUE)
        Aind[is.na(Aind)] <- 0L
        for (i in seq_len(nrow(na_ind))) {
            Aind[1, na_ind[i, 2]] <- Aind[1, na_ind[i, 2]] - 
                1L
        }
    }
    betahat <- quadprog::solve.QP.compact(as.matrix(Dmat), as.vector(dvec), 
        Amat, Aind)$sol
    RSS <- sum((y - (X %*% betahat))^2)
    if (trace) 
        print(paste("RSS:", RSS))
    attr(tree, "RSS") <- RSS
    bhat <- numeric(max(tree$edge))
    bhat[as.integer(lab)] <- betahat
    betahat <- bhat[tree$edge[, 1]] - bhat[tree$edge[, 2]]
    tree$edge.length <- betahat
    tree
}

nodeHeights = function (tree, ...) 
{
    if (hasArg(root.edge)) 
        root.edge <- list(...)$root.edge
    else root.edge <- FALSE
    if (root.edge) 
        ROOT <- if (!is.null(tree$root.edge)) 
            tree$root.edge
        else 0
    else ROOT <- 0
    nHeight <- function(tree) {
        tree <- reorder(tree)
        edge <- tree$edge
        el <- tree$edge.length
        res <- numeric(max(tree$edge))
        for (i in seq_len(nrow(edge))) res[edge[i, 2]] <- res[edge[i, 
            1]] + el[i]
        res
    }
    nh <- nHeight(tree)
    return(matrix(nh[tree$edge], ncol = 2L) + ROOT)
}

getTips = function(tree,node) {
  require(ape)
  if(node <= length(tree$tip.label)) {
    daughters <- tree$tip.label[node]
  } else {
    daughters <- extract.clade(tree, node = node)$tip.label
  }
  return(daughters)
}

#Note - only apply the "get_edge_from_tree" option if starting from an SNV only tree. Function will assume that all the existing edge length is SNVs.
correct_edge_length = function(node, tree, details, sensitivity_df, include_indels = TRUE, include_SNVs = TRUE, get_edge_from_tree=FALSE) {
  daughters <- getTips(tree = tree, node = node)
  #correct SNVs on edge, or set to 0 if want an indel only tree
  if(include_SNVs == TRUE) {
    if(get_edge_from_tree) {
    	nSNV=tree$edge.length[tree$edge[,2]==node]
    	} else {
    		nSNV = sum(details$node == node & details$Mut_type == "SNV")
    		}   		
    all_sens_SNVs <- sensitivity_df[sensitivity_df$Sample %in% daughters,"SNV_sensitivity"]
    branch_SNV_sens = 1 - prod(1-all_sens_SNVs)  
    new_nSNV = nSNV/branch_SNV_sens
  } else {
    new_nSNV <- 0
    }
  #correct INDELs on edge, or set to 0 if want an SNV only tree
  if(include_indels == TRUE) {
    nINDEL = sum(details$node == node & details$Mut_type == "INDEL")
    all_sens_INDELs <- sensitivity_df[sensitivity_df$Sample %in% daughters,"INDEL_sensitivity"]
    branch_INDEL_sens = 1 - prod(1-all_sens_INDELs)
    new_nINDEL = nINDEL/branch_INDEL_sens
  } else {
    new_nINDEL <- 0
  }
  new_edge_length = new_nSNV + new_nINDEL
  return(new_edge_length)
}

get_subset_tree = function(tree, details, v.field = "Mut_type", value = "SNV") {
  get_new_edge_length = function(node, tree, details,v.field,value) {
    sum(details$node == node & details[v.field] == value)
  }
  tree_subset = tree
  tree_subset$edge.length = sapply(tree$edge[,2], get_new_edge_length, tree = tree, details = details,v.field = v.field,value=value)
  return(tree_subset)
}

get_corrected_tree = function(tree, details, sensitivity_df, include_indels = TRUE, include_SNVs = TRUE,get_edge_from_tree=FALSE) {
  tree_c = tree
  tree_c$edge.length = sapply(tree$edge[,2], correct_edge_length, tree = tree, details = details, sensitivity_df = sensitivity_df, include_indels = include_indels, include_SNVs=include_SNVs,get_edge_from_tree=get_edge_from_tree)
  return(tree_c)
}

get_mut_burden = function(tree) {
  mut_burden = nodeHeights(tree)[tree$edge[,2] %in% 1:length(tree$tip.label),2]
}

get_mut_burden_stats = function(tree) {
  mut_burden = get_mut_burden(tree)
  cat(paste("Mean mutation burden is", round(mean(mut_burden),digits = 1),"\n"))
  cat(paste("Range of mutation burden is", round(range(mut_burden)[1],digits = 1),"to",round(range(mut_burden)[2],digits = 1),"\n"))
  cat(paste("Standard deviation of mutation burden is", round(sd(mut_burden),digits = 1),"\n"))
}
