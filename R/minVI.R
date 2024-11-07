### Finds optimal partition which minimizes the lower bound to the Variation of 
### Information obtained from Jensen's inequality where the expectation and 
### log are reversed.


minVI=function(psm, cls.draw=NULL, method=c("avg","comp","draws","all"), max.k=NULL){

	if(any(psm !=t(psm)) | any(psm >1) | any(psm < 0) | sum(diag(psm)) != nrow(psm) ){
	stop("psm must be a symmetric matrix with entries between 0 and 1 and 1's on the diagonals")}
    
	method <- match.arg(method, choices=method)
	if(method %in% c("draws","all") & is.null(cls.draw)) stop("cls.draw must be provided if method=''draws''")
	
	if(method == "avg" | method == "all"){
        if(is.null(max.k)) max.k <- ceiling(dim(psm)[1]/8)
        hclust.avg=stats::hclust(stats::as.dist(1-psm), method="average")
        cls.avg= t(apply(matrix(1:max.k),1,function(x) stats::cutree(hclust.avg,k=x)))
        VI.avg= VI.lb(cls.avg,psm)
        val.avg <- min(VI.avg)
        cl.avg <- cls.avg[which.min(VI.avg),] 
        if(method== "avg")  {
		output=list(cl=cl.avg, value=val.avg, method="avg")
		class(output)="c.estimate"
		return(output)
	  }
    }

    if(method == "comp" | method == "all"){
        if(is.null(max.k)) max.k <- ceiling(dim(psm)[1]/8)
        hclust.comp <- stats::hclust(stats::as.dist(1-psm), method="complete")
        cls.comp <-  t(apply(matrix(1:max.k),1,function(x) stats::cutree(hclust.comp,k=x)))
        VI.comp <- VI.lb(cls.comp,psm)
        val.comp <- min(VI.comp)
        cl.comp <- cls.comp[which.min(VI.comp),] 
        if(method== "comp")  {
		output=list(cl=cl.comp, value=val.comp, method="comp")
		class(output)="c.estimate"
		return(output)
	  }
    }

    if(method == "draws" | method == "all"){
		n=ncol(psm)
		EVI_lb_local=function(c){
			f=0
			for(i in 1:n){
				ind=(c==c[i])
				f=f+(log2(sum(ind))+log2(sum(psm[i,]))-2*log2(sum(ind*psm[i,])))/n
			}
			return(f)
		}
		VI.draws=apply(cls.draw,1,EVI_lb_local)
            val.draws <- min(VI.draws)
        	cl.draw <- cls.draw[which.min(VI.draws),] 
        	names(cl.draw) <- NULL
        	if(method== "draws") {
			output=list(cl=cl.draw, value=val.draws, method="draws")
			class(output)="c.estimate"
			return(output)
		}
    }

    vals <- c(val.avg, val.comp, val.draws)
    cls <- rbind(cl.avg,cl.comp,cl.draw)

    cls <- rbind(cls[which.min(vals),], cls)
    vals <- c(min(vals), vals)

    rownames(cls) <- names(vals) <- c("best","avg","comp","draws")
        colnames(cls) <- NULL    
        res <- list(cl=cls, value=vals, method="all")
	  
    class(res)="c.estimate"
    return(res)    
}
