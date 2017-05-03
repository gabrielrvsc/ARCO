#' Transforms a balanced panel into a list of matrices compatible with the fitArCo function 
#' 
#' Transforms a balanced panel into a list of matrices compatible with the fitArCo function. The user must identify the columns with the time, the unit identifier and the variables.
#' 
#' @param panel Balanced panel in a data.frame with columns for units and time.
#' @param time Name or index of the time column.
#' @param unit Name or index of the unit column.
#' @param variables Names or indexes of the columns containing the variables.
#' @export  
#' @examples 
#' # = Generate a small panel as example = #
#' set.seed(123)
#' time=sort(rep(1:100,2))
#' unit=rep(c("u1","u2"),100)
#' v1=rnorm(200)
#' v2=rnorm(200)
#' panel=data.frame(time=time,unit=unit,v1=v1,v2=v2)
#' head(panel)
#' 
#' data=panel_to_ArCo_list(panel,time="time",unit="unit",variables = c("v1","v2"))
#' head(data$v1)
#' 
#' @seealso \code{\link{fitArCo}}

panel_to_ArCo_list=function(panel,time,unit,variables){
  unit.nam=as.vector(unique(panel[,unit]))
  ArColist=list()
  for(q in 1:length(variables)){
    unitstore=matrix(rep(NA,nrow(panel)),ncol=length(unit.nam))
    colnames(unitstore)=unit.nam
    for(i in 1:length(unit.nam)){
      aux=panel[which(panel[,unit]==unit.nam[i]),c(time,unit,variables[q])]
      aux=aux[order(aux[,1]),]
      unitstore[,i]=aux[,3]
    }
    rownames(unitstore)=aux[,1]
    ArColist[[q]]=unitstore
  }
  names(ArColist)=variables
  return(ArColist)
}


