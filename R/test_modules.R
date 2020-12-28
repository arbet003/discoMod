###########################################
# similarity matrix for representing network structure -----------
###########################################
#'
#' @keywords internal
#'
#' @export

similarityMatrix=function(group1X,group2X,cortype="spearman",othertype=NULL,adjpower=1){
  #group1X: gene expression matrix for group 1, subjects in rows, genes in columns
  #group2X: gene expression matrix for group 2, subjects in rows, genes in columns
  #cortype="pearson","spearman", or "kendall"
  #othertype: "Adjacency" or "TOM": if othertype is not NULL, still need to specify the
  #cortype that will be converted to Adjacency or TOM
  #adjpower is the exponent on the adjacency matrix
  signed="signed" # if using othertype="Adjacency" or "TOM", only signed versions are considered.

  if(is.null(othertype)){
    similarityMat1=stats::cor(group1X,method=cortype)
    similarityMat2=stats::cor(group2X,method=cortype)
  } else {
    if(othertype=="Adjacency"){
      similarityMat1=WGCNA::adjacency(group1X,type=signed,power=adjpower)
      similarityMat2=WGCNA::adjacency(group2X,type=signed,power=adjpower)
    }
    if(othertype=="TOM"){
      adj1=WGCNA::adjacency(group1X,type=signed,power=adjpower)
      adj2=WGCNA::adjacency(group2X,type=signed,power=adjpower)

      similarityMat1=WGCNA::TOMsimilarity(adj1,TOMType = signed)
      similarityMat2=WGCNA::TOMsimilarity(adj2,TOMType = signed)
    }
  }

  colnames(similarityMat1)=colnames(group1X)
  rownames(similarityMat1)=colnames(group1X)
  colnames(similarityMat2)=colnames(group2X)
  rownames(similarityMat2)=colnames(group2X)

  ind=lower.tri(similarityMat1)
  V1=similarityMat1[ind]
  V2=similarityMat2[ind]

  ######stuff to return
  #similarityMat1: similarity matrix of modules for group 1
  #similarityMat2: similarity matrix of modules for group 2
  #V1: lower triangle of similarityMat1, used for test stats
  #V2  lower triangle of similarityMat2, used for test stats
  return(list(similarityMat1=similarityMat1,similarityMat2=similarityMat2,V1=V1,V2=V2,
              cortype=cortype,othertype=othertype,adjpower=adjpower,signed=signed))
}

########################################
# test statistics ----------------
########################################

#' @keywords internal
#'
#' @export
PND=function(V1,V2,p=6){
  diff= (abs(V1-V2))^p
  meandiff=mean(diff)
  teststat=meandiff^(1/p)
  return(teststat)
}

#' @keywords internal
#'
#' @export
GHD_test <- function(V1,V2) {
  GHD_diff <- mean(((V1-mean(V1))-(V2-mean(V2)))^2)
  teststats <- GHD_diff
  return(teststats)
}

#' @keywords internal
#'
#' @export
DispersionIndex=function(V1,V2){
  #V1: lower triangle of similarity matrix of module in group 1
  #V2: lower triangle of similarity matrix of module in group 2
  diff= (V1-V2)^2
  teststat=sqrt(mean(diff))
  return(teststat)
}

#' @keywords internal
#'
#' @export
MAD=function(V1,V2){
  diff= abs(V1-V2)
  teststat=mean(diff)
  return(teststat)
}

#' @keywords internal
#'
#' @export
pairedT=function(V1,V2){
  teststat=as.numeric(stats::t.test(V1,V2,paired=T)$stat)
  return(teststat)
}

#' @keywords internal
#'
#' @export
#'
wilcoxSRT=function(V1,V2){
  teststat=as.numeric(stats::wilcox.test(V1,V2,paired=T,exact=F)$stat)
  return(teststat)
}

######################################################################
# function to return all test stats ------------
######################################################################

#'
#' @keywords internal
#'
#' @export
#'
allTestStats=function(V1,V2){
  DIstat=DispersionIndex(V1,V2)
  MADstat=MAD(V1,V2)
  pairedTstat=pairedT(V1,V2)
  wilcoxSRTstat=wilcoxSRT(V1,V2)

  PND4stat=PND(V1,V2,p=4)
  PND6stat=PND(V1,V2,p=6)
  PND8stat=PND(V1,V2,p=8)
  PND20stat=PND(V1,V2,p=20)

  GHDstat <- GHD_test(V1,V2)

  allstats=c(PND4=PND4stat,
             PND6=PND6stat,
             PND8=PND8stat,
             PND20=PND20stat,
             GHD=GHDstat,
             DispIndex=DIstat,
             MAD=MADstat,
             pairedT=pairedTstat,
             wilcoxSRT=wilcoxSRTstat)
  return(allstats)
}

######################################################################
# permutations ----------
######################################################################

#' @keywords internal
#'
#' @export
#'
permTestStats=function(numperms,combinedGeneMat,similarityMatrix_object){
  #numperms: number of permutations
  #combinedGeneMat: (N1+N2)xP gene expression matrix,
  #     i.e. group 1 and group 2 are rbind into a single matrix
  #     needs column called "Group" that equals 1 for group 1 and 2 for group 2

  cortype=similarityMatrix_object$cortype
  adjpower=similarityMatrix_object$adjpower
  othertype=similarityMatrix_object$othertype

  permstats=vector()

  for(w in 1:numperms){
    permind=sample(1:nrow(combinedGeneMat),nrow(combinedGeneMat))
    permGeneMat=combinedGeneMat
    permGeneMat$Group=permGeneMat$Group[permind]

    permgroup1X=permGeneMat[which(permGeneMat$Group==1),]
    permgroup1X=permgroup1X[,-which(colnames(permgroup1X)=="Group")]

    permgroup2X=permGeneMat[which(permGeneMat$Group==2),]
    permgroup2X=permgroup2X[,-which(colnames(permgroup2X)=="Group")]

    permsimilarity=similarityMatrix(permgroup1X,permgroup2X,cortype=cortype,adjpower=adjpower,othertype=othertype)
    permV1=permsimilarity$V1 #lower triangle of similarity matrix in group 1
    permV2=permsimilarity$V2 #lower triangle of similarity matrix in group 2

    permallstats=allTestStats(permV1,permV2)
    permstats=rbind(permstats,permallstats)
  }
  permstats=data.frame(permstats)
  return(permstats)
}

#################################################
# permutation pvalues
#################################################

#' @keywords internal
#'
#' @export
#'
permpvalues=function(obs_stats,permstats){

  pvalues=vector()
  for(i in 1:length(obs_stats)){
    temp= permstats[,which(colnames(permstats)==names(obs_stats)[i])]

    # allow up to 5% missing values in permutations (very rarely occurs)
    checkNA= mean(is.na(temp))
    if(checkNA>0.05){
      pvalues[i] = NA
    } else {
      temp2= na.omit(temp)
      pvalues[i]=(sum(abs(temp2)>=abs(obs_stats[i]))+1)/(length(temp2)+1)
    }
  }
  names(pvalues)=names(obs_stats)
  return(pvalues)
}

#################################################
# wrapper functin for testing differentially coexpressed modules
#################################################

#' Test modules for differential co-expression
#'
#' For a given module, the null hypothesis is that the network structure is the same for
#' both phenotype groups, while the alternative hypothesis is that the network structure differs
#' between the two groups.
#' The following test statistics are calculated (all defined in a forthcoming manuscript),
#' using permutations to calculate p-values. PND: p-norm difference test (exponents 4,6,8, and 20),
#' GHD: Generalised Hamming Distance, DispIndex: Dispersion Index, MAD: mean absolute deviation,
#' pairedT: paired t-test, and wilcoxSRT: Wilcoxon signed rank test.
#' \strong{PND6} is the recommended default test based on simulation results from a forthcoming manuscript.
#'
#' @param group1_modules list of modules in phenotype group 1.  Each module (element of the list) is a matrix with samples in rows, and genes in columns.  You could use \code{find_modules()} to create this list, otherwise you can explore other clustering methods (e.g. hierarchical clustering) and create your own user-defined list
#' @param group2_modules list of modules in phenotype group 2.
#' @param cortype "spearman" (default), "pearson", or "kendall"; type of correlation used to measure the network structure within a module.
#' @param othertype If set to NULL (default), then correlation is used to represent the network.  "Adjacency" and "TOM" (topological overlap measure) are other options (only "signed" versions are considered) see \code{\link[WGCNA:adjacency]{WGCNA::adjacency()}} and \code{\link[WGCNA:TOMsimilarity]{WGCNA::TOMsimilarity()}} for details.  Option "TOM" is more computationally intensive.
#' @param adjpower if \code{othertype="Adjacency"} (or "TOM"), then \code{adjpower} represents the exponent used on all correlations (the higher the value, the more the correlations are shrunk towards zero).
#' @param numperms number of permutations used to calculate p-values
#' @param perm_stats TRUE or FALSE, whether or not to store all permutation test statistics (can take up a lot of memory if there are many modules and \code{numperms} is large)
#' @param parallel TRUE or FALSE, whether or not to create a doSNOW cluster used for parallel computing (the modules are evaluated in parallel, not the permutations)
#' @param cores number of cores to use if \code{parallel=T}.  Use \code{\link[parallel:detectCores]{parallel::detectCores()}} to check how many cores you have available.
#'
#' @return
#'
#' \itemize{
#'   \item \code{test_stats} : observed test statistics for each module (rows correspond to the order of modules in \code{group1_modules} and \code{group2_modules}).
#'   \item \code{pvalues} : p-values calculated using permutations.
#'   \item \code{qvalues} : false discovery rate multiple testing adjusted p-values (i.e. \code{stats::p.adjust()} with \code{method="fdr"})
#'   \item \code{perm_results} : if \code{perm_stats=TRUE} then \code{perm_results} will be a list of the permutation distribution of each test statistic for all modules.
#' }
#'
#' @import parallel
#' @import foreach
#' @import doSNOW
#'
#' @export
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

test_modules = function(group1_modules,
                        group2_modules,
                        cortype="spearman",
                        othertype=NULL,
                        adjpower=1,
                        numperms=500,
                        perm_stats=FALSE,
                        parallel=FALSE,
                        cores=2) {

  check = length(group1_modules)==length(group2_modules)
  if(check==F){
    stop("error: group1_modules and group2_modules must be lists of equal length, representing the same set of modules for both groups")
  }
  num_modules = sum(unlist(lapply(group1_modules, function(x) !is.null(x))))

  progress <- function(n, tag) if(tag%%cores==0){print(tag)}
  opts <- list(progress = progress)

  if(parallel==T){
    cl <- parallel::makeCluster(cores)
    doSNOW::registerDoSNOW(cl)
    parallel::clusterExport(cl,c(
      "allTestStats",
      "DispersionIndex",
      "GHD_test",
      "MAD",
      "pairedT",
      "permpvalues",
      "permTestStats",
      "PND",
      "similarityMatrix",
      "wilcoxSRT"
    ))
  } else {
    registerDoSEQ()
  }

  stats = foreach::foreach(w=1:num_modules,.options.snow = opts)%dopar%{
    start=Sys.time()

    group1X = group1_modules[[w]]
    group2X = group2_modules[[w]]

  if(!is.null(group1X) & !is.null(group2X)){
      ngenes_g1 = ncol(group1X)
      ngenes_g2 = ncol(group2X)

      if(ngenes_g1!=ngenes_g2){
        stop(paste0("error: the ",w,"th module has an unequal number of genes (columns) between the two groups"))
      }

      if(ngenes_g1>2){ # only test if >2 genes in module

          similarity=similarityMatrix(group1X, group2X,
                cortype=cortype,adjpower=adjpower,othertype=othertype)
          V1=similarity$V1 #lower triangle of similarity matrix in group 1
          V2=similarity$V2 #lower triangle of similarity matrix in group 2

          ###observed test statistics
          obs_stats=allTestStats(V1,V2)

          ###permutation test stats
          combinedGeneMat=rbind.data.frame(data.frame(group1X,Group=1),data.frame(group2X,Group=2))
          permstats=permTestStats(numperms = numperms,combinedGeneMat = combinedGeneMat,similarityMatrix_object = similarity)

          ###pvalues
          pvals=permpvalues(obs_stats,permstats)
          end=Sys.time()
          time=end-start

          if(perm_stats==F){
            permstats=NULL
          }

          results=list(test_stats=obs_stats,
                       pvalues=pvals,
                       permstats=permstats,
                       module=w,
                       time=time)
      } else {
        results=NULL
      }
  } else {
      results = NULL
  }
    print(w)
    return(results)
  }


  if(parallel==T){
    stopCluster(cl)
  }


  pvalues=as.data.frame(do.call(rbind,lapply(stats,function(x) x$pvalues)))

  qvalues=as.data.frame(apply(pvalues,2,function(x) stats::p.adjust(x,method="fdr")))

  test_stats=as.data.frame(do.call(rbind,lapply(stats,function(x) x$test_stats)))

  perm_results = lapply(stats,function(x) x$permstats)
  final_results = list(test_stats = test_stats,
                       pvalues = pvalues,
                       qvalues = qvalues,
                       perm_results=perm_results)
  return(final_results)
}

