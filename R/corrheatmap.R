#' Correlation heatmap for a single module
#'
#' Plot correlation heatmap of a module for both phenotype groups. \code{\link[ggcorrplot:ggcorrplot]{ggcorrplot::ggcorrplot()}} is used
#' to produce the heatmaps, which are returned in a list object in case the user
#' wants to further customize the plots using ggplot2.
#'
#' @param group1_mod matrix of a single module in phenotype group1
#' (samples in rows, genes in columns)
#' @param group2_mod matrix of the corresponding module in phenotype group1
#' (samples in rows, genes in columns)
#' @param cortype "spearman" (default), "pearson", or "kendall"; type of correlation
#' used to measure the network structure within a module.
#' @param hc.order TRUE/FALSE, whether to use hierarchical clustering to sort genes
#' @param group_order = 1 or 2, if hc.order=TRUE, then whether to perform the clustering
#' in group 1 or 2.  If you used \code{find_modules()} to derive the modules, then ideally,
#' \code{group_order} should be set to the same value that you set \code{cluster_group} too in \code{find_modules()}.
#' @param tl.cex font size for variable labels (if there are many variable names, e.g. > 100,
#' you might want to set \code{tl.cex = 0} to just remove the variable names from the plot)
#' @param ... additional parameters passed to \code{\link[ggcorrplot:ggcorrplot]{ggcorrplot::ggcorrplot()}}
#'
#' @return
#'
#' \itemize{
#'   \item \code{both_plots} : both heatmaps are combined into a single ggplot2 object using \code{ggpubr::ggarrange()}
#'   \item \code{group1_plot} : ggplot2 correlation heatmap of module in phenotype group 1
#'   \item \code{group2_plot} : ggplot2 correlation heatmap of module in phenotype group 2
#' }
#'
#' @export
#'
#'@examples
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

corrheatmap <- function(group1_mod,
                        group2_mod,
                        cortype = "spearman",
                        hc.order=TRUE,
                        group_order=1,
                        tl.cex=6, ...){
  check = all(colnames(group1_mod)==colnames(group2_mod))
  if(!check){
    stop("error: column names (genes) must be the same for both group1_mod and group2_mod")
  }


    if(hc.order==FALSE){
      g1_corr = stats::cor(group1_mod,method=cortype)
      g2_corr = stats::cor(group2_mod,method=cortype)
      p1 = ggcorrplot::ggcorrplot(g1_corr,hc.order=F,tl.cex=tl.cex,title="Group 1",...)
      p2 = ggcorrplot::ggcorrplot(g2_corr,hc.order=F,tl.cex=tl.cex,title="Group 2",...)
      figure = ggpubr::ggarrange(p1,p2,nrow=1,ncol=2,common.legend = T)
    } else {
      if(group_order==1){
        g1_corr = stats::cor(group1_mod,method=cortype)
        p1 = ggcorrplot::ggcorrplot(g1_corr,hc.order=T,tl.cex=tl.cex,title="Group 1",...)
        gorder = as.character(p1$data$Var1[1:ncol(group1_mod)])

        ind=vector()
        for(i in 1:length(gorder)){
          ind[i]=which(colnames(group2_mod)==gorder[i])
        }
        g2 = group2_mod[,ind]
        g2_corr = stats::cor(g2,method=cortype)
        p2 = ggcorrplot::ggcorrplot(g2_corr,hc.order=F,tl.cex=tl.cex,title="Group 2",...)
        figure = ggpubr::ggarrange(p1,p2,nrow=1,ncol=2,common.legend = T)
      }

      if(group_order==2){
        g2_corr = stats::cor(group2_mod,method=cortype)
        p2 = ggcorrplot::ggcorrplot(g2_corr,hc.order=T,tl.cex=tl.cex,title="Group 2",...)
        gorder = as.character(p2$data$Var1[1:ncol(group2_mod)])

        ind=vector()
        for(i in 1:length(gorder)){
          ind[i]=which(colnames(group1_mod)==gorder[i])
        }
        g1 = group1_mod[,ind]
        g1_corr = stats::cor(g1,method=cortype)
        p1 = ggcorrplot::ggcorrplot(g1_corr,hc.order=F,tl.cex=tl.cex,title="Group 1",...)
        figure = ggpubr::ggarrange(p1,p2,nrow=1,ncol=2,common.legend = T)
      }
    }
  return(list(both_plots = figure, group1_plot = p1, group2_plot = p2))
}
