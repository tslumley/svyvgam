## VGAM


svy_vglm<-function(formula,family,design,...){
    UseMethod("svy_vglm",design)
}

utils::globalVariables(".survey.prob.weights") 

svy_vglm.svyrep.design<-function(formula,family,design,...){
    
    vars<-intersect(all.vars(formula), colnames(design))
    surveydata<-model.frame(design)[,vars,drop=FALSE]
    
    pwts<-weights(design,"sampling")
    meanweight<-mean(pwts)
    surveydata$.survey.prob.weights<-pwts/meanweight
    surveydata$.survey.prob.weights[surveydata$.survey.prob.weights==0]<- 1e-9*max(pwts)

    
    fit0<-vglm(formula, family, data=surveydata, weights=.survey.prob.weights,...)

    theta<-coef(fit0)
    repwts<-weights(design, "analysis")
    thetas<-matrix(0,ncol=length(theta),nrow=ncol(repwts))
    for(i in 1:ncol(repwts)){
        surveydata$.survey.prob.weights<-repwts[,i]/meanweight
        surveydata$.survey.prob.weights[surveydata$.survey.prob.weights==0]<-1e-9*max(pwts)

        fit_i<-vglm(formula, family, data=surveydata, weights=.survey.prob.weights,...,coefstart=theta)
        thetas[i,]<-coef(fit_i)
    }

    v<-vcov(design,thetas, theta)
    dimnames(v)<-list(names(coef(fit0)),names(coef(fit0)))


    rval<-list(coef=coef(fit0), fit=fit0, design=design, var=v, call=sys.call())
    class(rval)<-c("svyrep_vglm","svy_vglm")
    rval
}


svy_vglm.survey.design<-function(formula, family, design,...){
    
    vars<-intersect(all.vars(formula), colnames(design))
    surveydata<-model.frame(design)[, vars, drop=FALSE]
    if(is.null(algorithm<-getOption("svyVGAM.algorithm")))
        algorithm<-1
    
    pwts<-weights(design,"sampling")
    meanweight<-mean(pwts)
    surveydata$.survey.prob.weights<-pwts
    surveydata$.survey.prob.weights[surveydata$.survey.prob.weights==0]<- 1e-9*max(pwts)
    
    fit<-vglm(formula, family, data=surveydata, weights=.survey.prob.weights,...)
    sfit<-summary(fit)
    
    naa<-fit@na.action
    if(!is.null(naa) && (length(naa)>0)){
        design<-design[-naa[[1]],]
        pwts<-pwts[-naa[[1]]]
    }
    if (algorithm == 1) {  # Original
        scores<-weights(fit, deriv = TRUE, type = "working")$deriv
        inv_inf<-sfit@cov.unscaled
        cons<-do.call(cbind, constraints(fit))
        mmat <- model.matrix(fit)
        case_index<-as.numeric(gsub("^(.+):.*", "\\1", rownames(mmat)))    
        mmatsum<- t(t(rowsum(mmat, case_index, reorder=FALSE ))/colSums(cons))
        inffuns<-(((scores/pwts)%*%cons)*mmatsum)%*% inv_inf
    } else { 
        ## based on advice from Thomas Yee
        ## however, fails for multinomial() family with only two categories
        dl.deta <- weights(fit, deriv = TRUE, type = "working")$deriv
        ## Remove prior weights
        dl.deta <- dl.deta / c(weights(fit, type = "prior"))  # use pwts?
        if (!is.matrix(dl.deta))
            dl.deta <- cbind(dl.deta)
        
        X.vlm <- model.matrix(fit, type = "vlm")
        nn <- nobs(fit)
        p.vlm <- ncol(X.vlm)
        M <- npred(fit)
        
        dl.dbeta.vlm <- matrix(0, nn, p.vlm)
        for (jay in 1:M) {
            vecTF <- rep(FALSE, M)
            vecTF[jay] <- TRUE  # Recycling
            dl.dbeta.vlm <- dl.dbeta.vlm +
                X.vlm[vecTF, , drop = FALSE] * dl.deta[, jay]
        }
        
        ## Any prior weights are in vcov(fit)
        inv_inf <- vcov(fit)  # Same as summary(fit)@cov.unscaled for VGLMs
        inffuns <- dl.dbeta.vlm %*% inv_inf
    }
    
    v<- vcov(svytotal(inffuns, design))
    dimnames(v)<-list(names(coef(fit)),names(coef(fit)))
    
    rval<-list(coef=coef(fit), fit=fit, var=v, naive.var=sfit@cov.unscaled*sfit@dispersion,
               design=design,algorithm=algorithm,call=sys.call())
    class(rval)<-"svy_vglm"
    rval
   
}

print.svy_vglm<-function(x, ...){
    print(x$design)
    show(x$fit)
    invisible(x)
}


vcov.svy_vglm<-function(object,...) object$var

coef.svy_vglm<-function(object,...) object$coef

halfp<- function(t) pnorm(-abs(t))/2
summary.svy_vglm<-function(object,...){
    object$coeftable<-cbind(Coef=coef(object), SE=SE(object), z=coef(object)/SE(object), p=halfp(coef(object)/SE(object)))
    class(object)<-"summary.svy_vglm"
    object
}

print.summary.svy_vglm<-function(x,...){
    print(x$call)
    print(x$design)
    printCoefmat(x$coeftable, has.Pvalue=TRUE, P.values=TRUE,signif.stars=FALSE)
    }
