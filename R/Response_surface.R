#' @import ggplot2
NULL

#' ResponseSurface.plot
#'
#' 
#'
#' @param rmap a ggplot2 plot object depicting a response map, produced by the function responseMap.
#' @param palette a character vector or two or more colors defining the colorscale to be applied to the fill variable in the plot.
#' @param cscale a value around which the colorscale applied to the fill variable will be centered. If NULL (the default), the colorscale will range from the minimum fill variable to the maximum.
#' @param xl The character string or expression used to label the x-axis
#' @param yl The character string or expression used to label the y-axis
#' @param zl The character string or expression used to label the colorscale representing the response variable. "Effect" by default.
#' @param title title for the plot
#' @return A ggplot2 plot object, with the appropriate scales and labels added.
#' @seealso braidReports::responseMap
#' @export
#' @examples 
#' \dontrun{
#' data(Dactolisib_Trametinib_combination)
#' growth_data<-net_growth_rate(Dactolisib_Trametinib_combination)
#' rmap <- braidReports::responseMap(Net_growth~CONC+CONC2,GD,logscale=T)
#'
#' ResponseSurface.plot(rmap,xl=expression(paste("[Dactolisib] (",mu,"M)", sep="")),
#'                     yl=expression(paste("[Trametinib] (",mu,"M)", sep="")),
#'                     zl="Net growth rate of \n sensitive cells",
#'                    palette=c("hotpink1","yellow","darkturquoise"))
#' }

ResponseSurface.plot <- function(rmap,palette=NULL,cscale=NULL,xl=expression(CONC),yl=expression(CONC2),zl="Effect",title="") {
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

#' DifferenceSurface.plot
#'
#' 
#'
#' @param rmap a ggplot2 plot object depicting a response map, produced by the function responseMap.
#' @param zcenter the value that forms the centerpoint (white) in the bidirectional color scale. If NULL (the default), the default setting for scale_fill_gradient2 is used.
#' @param xl The character string or expression used to label the x-axis
#' @param yl The character string or expression used to label the y-axis
#' @param zl The character string or expression used to label the colorscale representing the response variable. "Effect" by default.
#' @param low character corresponding to the color selected for the lowest values in the response map.
#' @param mid color for the mid point.
#' @param high character corresponding to the color selected for the highest values in the response map.
#' @param title title for the plot
#' @return A ggplot2 plot object, with the appropriate scales and labels added.
#' @seealso braidReports::responseMap
#' @export
#' @examples 
#' \dontrun{
#' data(Dactolisib_Trametinib_combination)
#' growth_data<-net_growth_rate(Dactolisib_Trametinib_combination)
#' GD=growth_data[,c('Cell.line','CONC','CONC2','Net_growth','Type','Birth_rate','Death_rate')]
#' GD=unique(GD)
#' rmap <- braidReports::responseMap(Net_growth~CONC+CONC2,GD,logscale=T,interpolate=FALSE)
#' DifferenceSurface.plot(rmap,zcenter=0.01,xl=expression(paste("[Dactolisib] (",mu,"M)", sep="")),
#' yl=expression(paste("[Trametinib] (",mu,"M)", sep="")),
#' zl="Net growth rate of \n sensitive cells",
#' mid="yellow",low="hotpink1",high="darkturquoise")
#'
#' }
DifferenceSurface.plot <- function(rmap,zcenter=NULL,xl=expression(CONC),yl=expression(CONC2),zl="Effect",low="#e9a3c9",mid="#f7f7f7",high="#4dac26",title=""){
  fp <- rmap+labs(x=xl,y=yl)+ggtitle(paste(title))
  if (!is.null(zcenter)) { fp <- fp+scale_fill_gradient2(zl,midpoint=zcenter,mid = mid,low=low,high=high) }
  #if (is.null(palette)) { palette=c("magenta","yellow","darkturquoise") }
  else { fp <- fp+scale_fill_gradient2(zl,mid = mid,low=low,high=high) }
  return(fp)
}
