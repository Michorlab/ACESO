#' @import ggplot2
NULL

#Response surface plot
#' @export
plot.ResponseSurface <- function(rmap,palette=NULL,cscale=NULL,xl=expression(CONC),yl=expression(CONC2),zl="Effect",title="") {
  fp <- rmap+labs(x=xl,y=yl)+ stat_contour(aes(z = rmap$data$z),colour="black", alpha=.2) +ggtitle(paste(title))
  if (is.null(palette)) { palette=c("magenta","yellow","darkturquoise") }
  if (!is.null(cscale)) {
    if (length(cscale)>1) { cscale <- cscale[1] }
    totmin <- min(rmap$data$z)
    totmax <- max(rmap$data$z)
    for (i in 1:length(rmap$layers)) {
      if (!is.null(rmap$layers[[i]]$data$z)) {
        totmin <- min(totmin,min(rmap$layers[[i]]$data$z))
        totmax <- max(totmax,max(rmap$layers[[i]]$data$z))
      }
    }
    dev <- max(totmax-cscale,cscale-totmin)+0.01
    clims <- c(cscale-dev,cscale+dev)
    fp <- fp + scale_fill_gradientn(zl,limits=clims,colours=palette)+theme_minimal()+ theme(text = element_text(size=14))
  } else { fp <- fp + scale_fill_gradientn(zl,colours=palette)+theme_minimal()+ theme(text = element_text(size=14)) }

  return(fp)
}

#' @export
plot.DifferenceSurface <- function(rmap,zcenter=NULL,xl=expression(CONC),yl=expression(CONC2),zl="Effect",low="#e9a3c9",mid="#f7f7f7",high="#4dac26") {
  fp <- rmap+labs(x=xl,y=yl)
  if (!is.null(zcenter)) { fp <- fp+scale_fill_gradient2(zl,midpoint=zcenter,mid = mid,low=low,high=high) }
  #if (is.null(palette)) { palette=c("magenta","yellow","darkturquoise") }
  else { fp <- fp+scale_fill_gradient2(zl,mid = mid,low=low,high=high) }
  return(fp)
}
