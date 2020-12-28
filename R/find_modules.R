
#' Find modules
#'
#' This is a wrapper function for \code{\link[mclust:Mclust]{mclust::Mclust()}}, which uses
#' Gaussian mixture models to find modules (clusters of correlated genes).
#' The Bayesian Information Criterion (BIC) is used to estimate the number of modules.
#' The modules are derived using only one of the phenotype groups, but the corresponding
#' modules are returned for both groups.  Later you can use the \code{test_modules()} function to
#' test whether the correlation structure within a module differs between the two groups.
#'
#' @param X1 matrix of genes for phenotype group 1 (rows are genes, columns are samples)
#' @param X2 matrix of genes for phenotype group 2 (rows are genes, columns are samples)
#' @param cluster_group which group you want to use for deriving the modules (= 1 or 2)
#' @param modelName the types of mclust covariance structures to consider.  See \code{\link[mclust:Mclust]{mclust::Mclust()}} for details.
#'   When clustering a large number of features, we recommend using diagonal covariance structures EII and VII.
#' @param num_modules the number of modules to consider.  Should be an integer vector specifying the number
#' of possible modules, e.g. 1:50, or you can specify an exact number of modules, e.g. num_modules=20.
#' @param ... additional arguments can be passed to \code{\link[mclust:Mclust]{mclust::Mclust()}}
#'
#' @return
#'
#' \itemize{
#'   \item \code{group1_modules} : list of modules in phenotype group 1 (from X1).  Each module is a matrix
#'   with samples in rows, and genes in columns.  Modules with < 3 genes are returned as NULL.
#'   \item \code{group2_modules} : list of modules in phenotype group 2 (from X2).
#'   \item \code{num_modules} : number of modules with 3 or more genes.
#'   \item \code{cluster_model} : the mclust model object for deriving the modules
#'   \item \code{module_classification} : how each row of X1,X2 is classified in terms of the resulting modules
#' }
#'
#'
#' @section References:
#'
#' Scrucca, L., Fop, M., Murphy, T. B., & Raftery, A. E. (2016). mclust 5: clustering, classification and
#' density estimation using Gaussian finite mixture models. The R journal, 8(1), 289.
#'
#'
#' @import mclust
#'
#'
#' @export
#'
#'
#' @examples
#'######## load example dataset
#'library(multtest)
#'data(golub)
#'X1 = golub[,which(golub.cl==0)]
#'X2 = golub[,which(golub.cl==1)]
#'rownames(X1) = golub.gnames[,3]
#'rownames(X2) = golub.gnames[,3]
#'
#'# use subset of 200 genes for example
#'set.seed(1234)
#'ind = sample(1:nrow(X1),200)
#'X1 = X1[ind,]
#'X2 = X2[ind,]
#'
#'######## Derive modules in group 1
#'modules = find_modules(X1,X2,cluster_group=1)
#'modules$num_modules # number of modules estimated by BIC (modules with < 3 genes are excluded)
#'ngm = unlist(lapply(modules$group1_modules, ncol)) # number of genes per module
#'summary(ngm)
#'
#'
#'######## test modules for differential co-expression
#'testmods = test_modules(group1_modules = modules$group1_modules,group2_modules = modules$group2_modules)
#'# View(testmods$pvalues)
#'# View(testmods$qvalues)
#'
#'which(testmods$pvalues$PND6 <= 0.05)
#'which(testmods$qvalues$PND6 <= 0.05)
#'
#'# use parallel computing:
#'# testmods = test_modules(group1_modules = modules$group1_modules,
#'#                        group2_modules = modules$group2_modules,
#'#                        parallel=TRUE,
#'#                        cores=4)
#'
#'
#'######## plot module 5
#'heat = corrheatmap(modules$group1_modules[[5]],modules$group2_modules[[5]])
#'# plot(heat$both_plots)
#'# plot(heat$group1_plot)
#'# plot(heat$group2_plot)

find_modules <- function(X1,X2,cluster_group=1, num_modules=1:100,
                         modelName = c("EII","VII"),...) {


    check = all(rownames(X1)==rownames(X2))
    if(check==F){
      stop("error: all rownames of X1 must equal rownames of X2, i.e. the rows for each group must represent the same genes")
    }

    if (cluster_group == "1") {
    cluster_model <- mclust::Mclust(X1, G=num_modules, modelName=modelName, ...)
    }

    if (cluster_group == "2") {
      cluster_model <- mclust::Mclust(X2, G=num_modules, modelName=modelName, ...)
    }

    moduleindex=cluster_model$classification
    group1_modules <- list()
    group2_modules <- list()
    for (i in 1:cluster_model$G) {
      ind=which(moduleindex==i)
      if(length(ind)>2){
        group1_modules[[i]]=t(X1[ind,])
        group2_modules[[i]]=t(X2[ind,])
      } else {
        group1_modules[[i]] = NULL
        group2_modules[[i]] = NULL
      }
    }
    # names(group1_modules)=paste0("module",1:cluster_model$G)
    # names(group2_modules)=paste0("module",1:cluster_model$G)

    num_modules = sum(unlist(lapply(group1_modules, function(x) !is.null(x))))

    results = list(group1_modules=group1_modules,
              group2_modules=group2_modules,
              num_modules = num_modules,
              cluster_model=cluster_model,
              module_classification=moduleindex)
    return(results)
}
