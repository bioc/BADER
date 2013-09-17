BADER <-
function (
    x,
    design,
    sizeFactors=TRUE,
    start=NULL,
    burn=1e3, reps=1e4,
    printEvery=1e2, saveEvery=1,
    t0 = 1e1,
    mode = "minimal" )
{
    #### check input ####
    if ( is.null(dim(x)))
        stop("Count table x should be a matrix")

    if ( !is.factor(design))
        stop("'design' parameter should be a factor")

    if ( dim(x)[2] > dim(x)[1])
        warning("The count table x has more rows than columns. In a typical
            usage example where columns represent samples and rows features
            the contrary should be the case")

    if ( !is.null(start) & typeof(start) != "list")
        stop("The 'start' parameter needs to be either NULL or a list of starting parameters")

    if ( length(levels(design)) != 2 )
        stop(paste("The 'design' parameter contains",length(levels(design)),
            "groups. At the moment, BADER can only handle two treatment groups"))
            
    if ( length(design) != dim(x)[2] )
        stop("The length of the 'design' parameter does not match the number
            of samples in the count table x")

    match.arg(mode,c("minimal","full","extended"))
    ##


    #### initialize helper variables ####
    groups <- levels(design)
    
    # indices of two groups
    iA <- design == groups[1]
    iB <- design == groups[2]

    nA <- sum(iA)
    nB <- sum(iB)
    m <- dim(x)[1]
   
    #### estimate size factors ####
    if ( sizeFactors == FALSE )
    {
        sA <- rep(1, nA)
        sB <- rep(1, nB)
    }

    if ( sizeFactors == TRUE )
    {
        geoMeans <- exp(rowMeans(log(x)))
        sfac <- apply(x,2,function(x)
            median( (x / geoMeans)[geoMeans>0]))

        sA <- sfac[iA]
        sB <- sfac[iB]
    }

    #### guess starting parameters for MCMC ####    
    ## use simple method-of-moments estimators ##  
    kA <- t(x[,iA])
    kB <- t(x[,iB])

    lambdaA <- ifelse(kA/sA > 0, log(kA/sA), NA)
    lambdaA[is.na(lambdaA)] <- min(lambdaA, na.rm=TRUE)

    lambdaB <- ifelse(kB/sB > 0, log(kB/sB), NA)
    lambdaB[is.na(lambdaB)] <- min(lambdaB, na.rm=TRUE)

    muA <- colMeans(lambdaA)
    muB <- colMeans(lambdaB)

    gam <- muB - muA

    pi0 <- 0.1
    ind <- rep(0,m)
    cutoff <- sort(abs(gam))[round((1-pi0)*m)] 
    ind[abs(gam) > cutoff] <- 1
    gam[ind == 0] <- 0
    sigmaGamma <- var(gam[gam!=0])
  
    varA <- apply(kA/sA,2,var)
    varB <- apply(kB/sB,2,var)

    alphaA <- rep(NA,m)
    alphaB <- rep(NA,m)
    alphaA[varA-exp(muA) > 0] <- log(((varA - exp(muA))*exp(-2*muA))[varA-exp(muA) > 0])
    alphaB[varB-exp(muB) > 0] <- log(((varB - exp(muB))*exp(-2*muB))[varB-exp(muB) > 0])

    tau <- var(c(alphaA,alphaB), na.rm=TRUE)
    alphaA[is.na(alphaA)] <- mean(alphaA, na.rm=TRUE)
    alphaB[is.na(alphaB)] <- mean(alphaB, na.rm=TRUE)


    psi0 <- mean(c(alphaA,alphaB))
    
    #### use start parameters if present ####
    if ( typeof(start) == "list")
    {
        params <- c("lambdaA", "lambdaB", "muA", "gam", "ind", 
            "alphaA", "alphaB", "psi0", "tau", "sigmaGamma", "pi0")

        for ( param in intersect(names(start),params) )
        {
            if ( !is.null(dim(get(param))) && dim(get(param)) != dim(start[[param]]))
                stop(paste("Starting parameter '",param,"' does not have the correct dimensions",sep=""))
            if ( length(get(param)) != length(start[[param]]))
                stop(paste("Starting parameter '",param,"' does not have the correct length",sep=""))
            
            assign(param,start[[param]])        
        }
    }

    #### perform MCMC ####
    if ( mode == "minimal" )
    {
        results <- .C("rnaseq", as.double(kA), as.double(kB), as.double(sA), as.double(sB),
            as.integer(nA), as.integer(nB), as.integer(m), as.integer(burn), as.integer(reps),
            as.integer(saveEvery), as.integer(printEvery), as.integer(t0),
            as.double(lambdaA), as.double(lambdaB), as.double(ind), as.double(muA),
            as.double(gam), as.double(alphaA), as.double(alphaB), as.double(pi0),
            as.double(sigmaGamma), as.double(psi0), as.double(tau), as.integer(0)
        )
    }
    if ( mode == "full" )
    {
        distGamma <- rep(0,floor(reps/saveEvery)*m)
        results <- .C("rnaseq_post_dist", as.double(kA), as.double(kB), as.double(sA), as.double(sB),
            as.integer(nA), as.integer(nB), as.integer(m), as.integer(burn), as.integer(reps),
            as.integer(saveEvery), as.integer(printEvery), as.integer(t0),
            as.double(lambdaA), as.double(lambdaB), as.double(ind), as.double(muA),
            as.double(gam), as.double(alphaA), as.double(alphaB), as.double(pi0),
            as.double(sigmaGamma), as.double(psi0), as.double(tau),
            as.double(distGamma), as.integer(0)
        )
    }
    if ( mode == "extended" )
    {
        distGamma <- distMuA <- distAlphaA <- distAlphaB <- rep(0,floor(reps/saveEvery)*m)
        distPi0 <- distSigmaGamma <- distPsi0 <- distTau <- rep(0,floor(reps/saveEvery))
        
        results <- .C("rnaseq_verbose", as.double(kA), as.double(kB), as.double(sA), as.double(sB),
            as.integer(nA), as.integer(nB), as.integer(m), as.integer(burn), as.integer(reps),
            as.integer(saveEvery), as.integer(printEvery), as.integer(t0),
            as.double(lambdaA), as.double(lambdaB), as.double(ind), as.double(muA),
            as.double(gam), as.double(alphaA), as.double(alphaB), as.double(pi0),
            as.double(sigmaGamma), as.double(psi0), as.double(tau),
            as.double(distGamma), as.double(distMuA), as.double(distAlphaA), as.double(distAlphaB),
            as.integer(0)
        )
    }

    #### return results ####
    if ( mode == "minimal" )
    {
        if ( results[[24]] == 0 )
        {
            ### sampling was not successfull
            warning("Sampling was not sucessful. Maybe there are zero columns in the data? Also,
                  changing the t0 parameter could help")
        }
        obj <- list (
            logMeanA = results[[16]],
            logMeanB = (results[[16]] + results[[17]]),
            logFoldChange = results[[17]],
            logDispA = results[[18]],
            logDispB = results[[19]],
            diffProb = results[[15]]
        )
        class (obj) <- c(class(obj),"BADER") 
    }
    if ( mode == "full" )
    {
        if ( results[[25]] == 0 )
        {
            ### sampling was not successfull
            warning("Sampling was not sucessful. Maybe there are zero columns in the data? Also,
                  changing the t0 parameter could help")
        }

        obj <- list (
            logMeanA = results[[16]],
            logMeanB = (results[[16]] + results[[17]]),
            logFoldChange = results[[17]],
            logDispA = results[[18]],
            logDispB = results[[19]],
            diffProb = results[[15]],
            logFoldChangeDist = matrix(results[[24]],ncol=m)
        )
        class (obj) <- c(class(obj),"BADER") 
    }
    if ( mode == "extended" )
    {
        if ( results[[28]] == 0 )
        {
            ### sampling was not successfull
            warning("Sampling was not sucessful. Maybe there are zero columns in the data? Also,
                  changing the t0 parameter could help")
        }

        obj <- list (
            logMeanA = results[[16]],
            logMeanB = (results[[16]] + results[[17]]),
            logFoldChange = results[[17]],
            logDispA = results[[18]],
            logDispB = results[[19]],
            diffProb = results[[15]],
            logFoldChangeDist = matrix(results[[24]],ncol=m),
            logMeanADist = matrix(results[[25]],ncol=m),
            logMeanBDist = matrix(results[[24]] + results[[25]],ncol=m),
            logDispADist = matrix(results[[26]],ncol=m),
            logDispBDist = matrix(results[[27]],ncol=m)
        )
        class (obj) <- c(class(obj),"BADER") 
    }

    return(obj)
}
