#############################################################################
### Identify doublets in single-cell RNA-seq data, using Scrublet method. ###
###              Chengxiang Qiu, Genome Sciences, UW, 2018-11-18          ###
#############################################################################

###################################################
### The main function used to identify doublets ###
###################################################

#' The main function used to identify doublets
#' @param exp expression matrix of a single-cell RNA-seq data set (gene, cell)
#' @param n_neighbors the number of neighbors used to create the KNN graph
#' @param doublet_score_threshold the threshold used to separate doublets
#' @param sim_doublet_ratio the expected doublet ratio
#' @param stdev_doublet_rate the standard deviation of doublet rate
#' @param synthetic_doublet_umi_subsampling the subsampling rate while simulation
#' @param use_approx_neighbors whether use approx neighbors
#' @param distance_metric which method use measure distance of two nodes
#' @param min_counts minimum read count of cell
#' @param min_cells minimum detected cell of gene
#' @param min_gene_variability_pctl minumum gene variability percentage
#' @param log_transform whether perform log transformation
#' @param z_score whether transfer to z_score
#' @param n_prin_comps the number of prin comps used
#' @param verbose whether display processing message
#' @return the doublet label and score for each cell
#' @export
scrubDoublets <- function(exp,
                        n_neighbors = NULL,
                        doublet_score_threshold = NULL,
                        sim_doublet_ratio = 2.0,
                        expected_doublet_rate = 0.06,
                        stdev_doublet_rate = 0.02,
                        synthetic_doublet_umi_subsampling = 1.0,
                        use_approx_neighbors = TRUE,
                        distance_metric = 'euclidean',
                        min_counts = 2,
                        min_cells = 3,
                        min_gene_variability_pctl = 85,
                        log_transform = FALSE,
                        z_score = TRUE,
                        n_prin_comps = 30,
                        verbose = TRUE){

    if (is.null(n_neighbors)) n_neighbors <- round(0.5 * sqrt(ncol(exp)))

    if (verbose) message("Preprocessing...")
    E_obs <- t(exp) ### E_obs, ncell * ngene
    total_counts_obs <- apply(E_obs, 1, sum)

    E_obs_norm <- pipeline_normalize(E_obs, total_counts_obs)
    gene_filter <- pipeline_get_filter(E_obs_norm)
    E_obs <- E_obs[,gene_filter]
    E_obs_norm <- E_obs_norm[,gene_filter]

    if (verbose) message("Simulating doublets...")
    simulateDoublets.res <- simulateDoublets(E_obs, total_counts_obs, sim_doublet_ratio, synthetic_doublet_umi_subsampling)
    E_sim <- simulateDoublets.res$E_sim
    total_counts_sim <- simulateDoublets.res$total_counts_sim
    E_obs_norm <- pipeline_normalize(E_obs, total_counts_obs, postnorm_total = 1e6)
    E_sim_norm <- pipeline_normalize(E_sim, total_counts_sim, postnorm_total = 1e6)

    if (log_transform) {
        E_obs_norm <- pipeline_log_transform(E_obs_norm)
        E_sim_norm <- pipeline_log_transform(E_sim_norm)
    }

    if (z_score) {
        gene_mean <- apply(E_obs_norm, 2, mean)
        gene_std <- apply(E_obs_norm, 2, sd)
        E_obs_norm <- pipeline_zscore(E_obs_norm, gene_mean, gene_std)
        E_sim_norm <- pipeline_zscore(E_sim_norm, gene_mean, gene_std)
    }

    pca.res <- pipeline_pca(E_obs_norm, E_sim_norm, n_prin_comps)

    if (verbose) message("Calculating doublet scores...")
    doublet_scores <- calculateDoubletScores(pca.res$pca_obs, pca.res$pca_sim, n_neighbors)

    if (is.null(doublet_score_threshold)) {
        if (verbose) message("Histogram of doublet scores...")
        predicted_threshold <- histogramDoubletScores(doublet_scores$doublet_scores_obs, doublet_scores$doublet_scores_sim)
        doublet_score_threshold <- predicted_threshold
    }

    if (verbose) message("Call transcriptomes as doublets...")
    predicted_doublets <- callDoublets(doublet_scores$doublet_scores_obs, doublet_scores$doublet_scores_sim, expected_doublet_rate, doublet_score_threshold, verbose)

    return(list(scrubDoublets = predicted_doublets, doublet_scores_obs = doublet_scores$doublet_scores_obs, doublet_scores_sim = doublet_scores$doublet_scores_sim))

}

##################################################
### manually reset the doublet_score_threshold ###
##################################################

scrubDoublets_resetThreshold <- function(scrubDoublets_res,
                                doublet_score_threshold = NULL,
                                verbose = TRUE){

    if(is.null(doublet_score_threshold)) message("Please set doublet_score_threshold.")
    else{
        Ld_obs <- scrubDoublets_res$doublet_scores_obs
        Ld_sim <- scrubDoublets_res$doublet_scores_sim
        threshold <- doublet_score_threshold

        predicted_doublets <- Ld_obs > threshold

        detected_doublet_rate <- sum(Ld_obs > threshold) / length(Ld_obs)
        detectable_doublet_fraction <- sum(Ld_sim>threshold) / length(Ld_sim)
        overall_doublet_rate <- detected_doublet_rate / detectable_doublet_fraction

        if (verbose) {
            message(paste0("Set threshold at doublet score = ", round(doublet_score_threshold,2)))
            message(paste0("Detected doublet rate = ", 100*round(detected_doublet_rate,4), "%"))
            message(paste0("Estimated detectable doublet fraction = ", 100*round(detectable_doublet_fraction,4), "%"))
            message("Overall doublet rate:")
            message(paste0("Estimated  = ", 100*round(overall_doublet_rate,4), "%"))
        }

        return(list(scrubDoublets = predicted_doublets, doublet_scores_obs = Ld_obs, doublet_scores_sim = Ld_sim))

    }

}


############################
### pipeline: preprocess ###
############################

pipeline_normalize <- function(E, total_counts, postnorm_total = NULL){

    if (is.null(postnorm_total)){
        total_counts_mean <- mean(total_counts)
    } else {
        total_counts_mean <- postnorm_total
    }

    ncell <- nrow(E)
    w <- matrix(0, ncell, ncell)
    diag(w) <- total_counts_mean / total_counts
    Enorm <- w %*% E

    return(Enorm)
}

pipeline_get_filter <- function(Enorm, min_counts = 3, min_cells = 3, min_gene_variability_pctl = 85){

    vscores.res <- get_vscores(Enorm)
    ix2 <- vscores.res$Vscores > 0
    Vscores <- vscores.res$Vscores[ix2]
    gene_ix <- vscores.res$gene_ix[ix2]
    mu_gene <- vscores.res$mu_gene[ix2]
    FF_gene <- vscores.res$FF_gene[ix2]
    min_vscore <- quantile(Vscores, min_gene_variability_pctl/100)
    ix <- (apply(Enorm[,gene_ix]>=min_counts, 2, sum) >= min_cells) & (Vscores >= min_vscore)

    return(gene_ix[ix])
}

#' @importFrom neldermead fminsearch
get_vscores <- function(Enorm, min_mean = 0, nBins = 50, fit_percentile = 0.1, error_wt = 1){

    ncell <- nrow(Enorm)
    mu_gene <- apply(Enorm, 2, mean)
    gene_ix <- c(1:ncol(Enorm))[mu_gene > min_mean]
    mu_gene <- mu_gene[gene_ix]

    tmp <- Enorm[,gene_ix]
    tmp <- tmp^2
    var_gene <- apply(tmp, 2, mean) - mu_gene^2
    FF_gene <- var_gene / mu_gene

    data_x <- log(mu_gene)
    data_y <- log(FF_gene / mu_gene)

    tmp <- runningquantile(data_x, data_y, fit_percentile, nBins)
    x <- tmp$xOut[!is.na(tmp$yOut)]
    y <- tmp$yOut[!is.na(tmp$yOut)]

    gLog <- function(x0, x1, x2) log(x1 * exp(-x0) + x2)
    tmp <- log(FF_gene[mu_gene>0])
    tmp <- hist(tmp, breaks=seq(min(tmp), max(tmp), l=201), plot=F)
    h <- tmp$counts
    b <- tmp$breaks
    b <- b[-201] + diff(b)/2
    max_ix <- which.max(h)
    c <- max(exp(b[max_ix]), 1)
    errFun <- function(b2) sum(abs(gLog(x, c, b2)-y)^error_wt)
    b0 <- 0.1
    b <- neldermead::fminsearch(errFun, b0)$simplexopt$x[1,]
    a <- c / (1+b) - 1

    v_scores <- FF_gene / ((1+a)*(1+b) + b*mu_gene)

    return(list(Vscores=v_scores, gene_ix=gene_ix, mu_gene=mu_gene, FF_gene=FF_gene, a=a, b=b))
}

runningquantile <- function(x, y, p, nBins){

    ind <- order(x)
    x <- x[ind]
    y <- y[ind]

    dx <- (x[length(x)] - x[1]) / nBins
    xOut <- seq(x[1]+dx/2, x[length(x)]-dx/2, len=nBins)
    yOut <- rep(0, length(xOut))

    for (i in 1:length(xOut)){
        ind <- (x >= (xOut[i]-dx/2)) & (x < (xOut[i]+dx/2))
        if (sum(ind)>0){
            yOut[i] <- quantile(y[ind], p/100)
        }
        else{
            if (i>1){
                yOut[i] <- yOut[i-1]
            }
            else {
                yOut[i] <- NA
            }
        }
    }

    return(list(xOut=xOut, yOut=yOut))

}

pipeline_log_transform <- function(E, pseudocount = 1){
    X <- log10(E + pseudocount)
    return(X)
}

pipeline_zscore <- function(E, gene_mean, gene_std){

    E <- t(sweep(E, 2, gene_mean))

    nrow <- nrow(E)
    w <- matrix(0, nrow, nrow)
    diag(w) <- 1 / gene_std
    X <- w %*% E

    return(t(X))
}

pipeline_pca <- function(X_obs, X_sim, n_prin_comps){

    pca <- prcomp(X_obs, rank.=n_prin_comps)
    pca_obs <- pca$x
    pca_sim <- scale(X_sim, pca$center, pca$scale) %*% pca$rotation

    return(list(pca_obs = pca_obs, pca_sim = pca_sim))
}

#################################################
### Simulate doublets from observed read cout ###
#################################################

simulateDoublets <- function(E_obs,
                            tatal_counts_obs,
                            sim_doublet_ratio = 2.0,
                            synthetic_doublet_umi_subsampling = 1.0){

    n_obs <- nrow(E_obs)
    n_sim <- round(n_obs * sim_doublet_ratio)

    pair_ix <- matrix(,n_sim,2)
    for(i in 1:n_sim){
        pair_ix[i,] <- sample(1:n_obs,2)
    }

    E1 <- E_obs[pair_ix[,1],]
    E2 <- E_obs[pair_ix[,2],]
    tots1 <- tatal_counts_obs[pair_ix[,1]]
    tots2 <- tatal_counts_obs[pair_ix[,2]]

    if (synthetic_doublet_umi_subsampling < 1){
        simulateDoublets.tmp <- subsampleCounts(E1 + E2, synthetic_doublet_umi_subsampling, tots1+tots2)
        E_sim <- simulateDoublets.tmp[[1]]
        total_counts_sim <- simulateDoublets.tmp[[2]]
    } else {
        E_sim <- E1 + E2
        total_counts_sim <- tots1 + tots2
    }

    return(list(E_sim = E_sim, total_counts_sim = total_counts_sim, pair_ix = pair_ix))

}

subsampleCounts <- function(E, rate, original_totals){
    E <- matrix(rbinom(nrow(E) * ncol(E),round(E),rate),nrow(E),ncol(E))
    current_totals <- apply(E, 1, sum)
    unsampled_orig_totals <- original_totals - current_totals
    unsampled_downsamp_totals <- rbinom(length(unsampled_orig_totals), round(unsampled_orig_totals), rate)
    final_downsamp_totals <- current_totals + unsampled_downsamp_totals

    return(list(E, final_downsamp_totals))
}

################################
### Calculate doublet scores ###
################################
#' @importFrom RANN nn2
calculateDoubletScores <- function(pca_obs,
                                pca_sim,
                                n_neighbors,
                                expected_doublet_rate = 0.06,
                                stdev_doublet_rate = 0.02,
                                distance_metric = "euclidean"){

    n_obs <- nrow(pca_obs)
    n_sim <- nrow(pca_sim)
    manifold <- rbind(pca_obs, pca_sim)
    doub_labels <- c(rep(0, n_obs), rep(1, n_sim))

    # Find k_adj nearest neighbors
    k_adj <- round(n_neighbors * (1+n_sim/n_obs))
    #if (distance_metric %in% c("euclidean")) neighbors <- get.knn(manifold, k = k_adj)$nn.index
    if (distance_metric %in% c("euclidean")) neighbors <- RANN::nn2(manifold,k = k_adj)$nn.idx[,-1]

    # Calculate doublet score based on ratio of simulated cell neighbors vs. observed cell neighbors
    doub_neigh_mask <- matrix(doub_labels[neighbors] == 1, nrow(neighbors), ncol(neighbors))
    n_sim_neigh <- apply(doub_neigh_mask, 1, sum)
    n_obs_neigh <- k_adj - n_sim_neigh

    rho <- expected_doublet_rate
    r <- n_sim / n_obs
    nd <- n_sim_neigh
    ns <- n_obs_neigh
    N <- k_adj

    # Bayesian
    q <- (nd+1)/(N+2)
    Ld <- q*rho/r/(1-rho-q*(1-rho-rho/r))

    se_q <- sqrt(q*(1-q)/(N+3))
    se_rho <- stdev_doublet_rate

    se_Ld <- q*rho/r / (1-rho-q*(1-rho-rho/r))**2 * sqrt((se_q/q*(1-rho))**2 + (se_rho/rho*(1-q))**2)

    doublet_scores_obs <- Ld[doub_labels == 0]
    doublet_scores_sim <- Ld[doub_labels == 1]
    doublet_errors_obs <- se_Ld[doub_labels==0]
    doublet_errors_sim <- se_Ld[doub_labels==1]

    return(list(doublet_scores_obs = doublet_scores_obs, doublet_scores_sim = doublet_scores_sim, doublet_errors_obs = doublet_errors_obs, doublet_errors_sim = doublet_errors_sim))

}

###################################################################
### Plot the histograme for doublet scores and detect threshold ###
###################################################################

#' @import ggplot2
#' @importFrom gridExtra grid.arrange
histogramDoubletScores <- function(doublet_scores_obs, doublet_scores_sim){

    # estimte the threshold based on kmeans cluster
    km <- kmeans(doublet_scores_sim, centers=2)
    clust <- as.factor(km$cluster)
    predicted_threshold <- (max(doublet_scores_sim[clust==1]) + min(doublet_scores_sim[clust==2]))/2

    dat_obs <- data.frame(doublet_scores = doublet_scores_obs, clust = rep(1, length(doublet_scores_obs)))
    dat_obs$clust[dat_obs$doublet_scores > predicted_threshold] <- 2
    dat_obs$clust <- factor(dat_obs$clust)

    dat_sim <- data.frame(doublet_scores = doublet_scores_sim, clust = rep(1, length(doublet_scores_sim)))
    dat_sim$clust[dat_sim$doublet_scores > predicted_threshold] <- 2
    dat_sim$clust <- factor(dat_sim$clust)

    p_obs <- ggplot(dat_obs, aes(x = doublet_scores))
    p_obs <- p_obs + geom_histogram(aes(fill = clust), binwidth = 0.02, color = "grey50")
    p_obs <- p_obs + geom_vline(xintercept = predicted_threshold, color = "blue")
    p_obs <- p_obs + labs(x="Doublet scores", y="Counts", title="Observed Cells") + theme_classic(base_size = 10) + theme(legend.position="none") + theme(plot.title = element_text(hjust = 0.5))

    p_sim <- ggplot(dat_sim, aes(x = doublet_scores))
    p_sim <- p_sim + geom_histogram(aes(fill = clust), binwidth = 0.02, color = "grey50")
    p_sim <- p_sim + geom_vline(xintercept = predicted_threshold, color = "blue")
    p_sim <- p_sim + labs(x="Doublet scores", y="Counts", title="Simulated Doublets") + theme_classic(base_size = 10) + theme(legend.position="none") + theme(plot.title = element_text(hjust = 0.5))

    p_obs2 <- ggplot(dat_obs, aes(x = doublet_scores)) + stat_density(geom="line", color="red") + geom_vline(xintercept = predicted_threshold, color = "blue")
    p_obs2 <- p_obs2 + labs(x="Doublet scores", y="Density", title="") + theme_classic(base_size = 10) + theme(legend.position="none") + theme(plot.title = element_text(hjust = 0.5))
    p_sim2 <- ggplot(dat_sim, aes(x = doublet_scores)) + stat_density(geom="line", color="red") + geom_vline(xintercept = predicted_threshold, color = "blue")
    p_sim2 <- p_sim2 + labs(x="Doublet scores", y="Density", title="") + theme_classic(base_size = 10) + theme(legend.position="none") + theme(plot.title = element_text(hjust = 0.5))

    pdf("histogram of doublet scores.pdf",8,8)
    gridExtra::grid.arrange(p_obs, p_sim, p_obs2, p_sim2, nrow = 2, ncol = 2)
    dev.off()

    return(predicted_threshold)

}



###################################################
### Call transcriptomes as doublets or singlets ###
###################################################

#' reset the doublet score threshold
#' @param doublet_scores_obs doublet scores of obs cells
#' @param doublet_scores_sim doublet scores of sim cells
#' @param expected_doublet_rate expected doublet rate
#' @param doublet_score_threshold reset the threshold of doublet score
#' @param verbose display processing messages
#' @return predicted doublet scores
#' @export
callDoublets <- function(doublet_scores_obs,
                        doublet_scores_sim,
                        expected_doublet_rate,
                        doublet_score_threshold,
                        verbose){

    Ld_obs <- doublet_scores_obs
    Ld_sim <- doublet_scores_sim
    threshold <- doublet_score_threshold

    predicted_doublets <- Ld_obs > threshold

    detected_doublet_rate <- sum(Ld_obs > threshold) / length(Ld_obs)
    detectable_doublet_fraction <- sum(Ld_sim>threshold) / length(Ld_sim)
    overall_doublet_rate <- detected_doublet_rate / detectable_doublet_fraction

    if (verbose) {
        message(paste0("Set threshold at doublet score = ", round(doublet_score_threshold,2)))
        message(paste0("Detected doublet rate = ", 100*round(detected_doublet_rate,4), "%"))
        message(paste0("Estimated detectable doublet fraction = ", 100*round(detectable_doublet_fraction,4), "%"))
        message("Overall doublet rate:")
        message(paste0("Expected   = ", 100*round(expected_doublet_rate,4), "%"))
        message(paste0("Estimated  = ", 100*round(overall_doublet_rate,4), "%"))
    }

    return(predicted_doublets)

}

