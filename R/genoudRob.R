#
#  multinomRob
#
#  Walter R. Mebane, Jr.
#  Cornell University
#  http://macht.arts.cornell.edu/wrm1/
#  wrm1@macht.arts.cornell.edu
#
#  Jasjeet Singh Sekhon 
#  Harvard University
#  http://jsekhon.fas.harvard.edu/
#  jsekhon@fas.harvard.edu
#
#  $Id: genoudRob.R,v 1.3 2004/02/14 22:21:35 wrm1 Exp $
#
###################################
#New Front End for Genoud, with tuned defaults
###################################

#sets genoud.parms defaults
genoudParms  <- function(genoud.parms)
  {
    #set user controlled defaults
    if (is.null(genoud.parms$pop.size))
      genoud.parms$pop.size  <- 1000;

    if (is.null(genoud.parms$max.generations))
      genoud.parms$max.generations  <- 100;
    
    if (is.null(genoud.parms$wait.generations))
      genoud.parms$wait.generations  <- 10;
    
    if (is.null(genoud.parms$hard.generation.limit))
      genoud.parms$hard.generation.limit  <- FALSE;

    #this is redundant, but maintains clarity  
    if (is.null(genoud.parms$MemoryMatrix))
      genoud.parms$MemoryMatrix  <- NULL;
  
    if (is.null(genoud.parms$Debug))
      genoud.parms$Debug  <- FALSE ;

    #this is redundant, but maintains clarity   
    if (is.null(genoud.parms$Domains))
      genoud.parms$Domains  <- NULL;

    if (is.null(genoud.parms$scale.domains))
      genoud.parms$scale.domains  <- 10;
    
    if (is.null(genoud.parms$boundary.enforcement))
      genoud.parms$boundary.enforcement  <- 0;
    
    if (is.null(genoud.parms$solution.tolerance))
      genoud.parms$solution.tolerance  <- 0.0000001;
    
    if (is.null(genoud.parms$BFGS))
      genoud.parms$BFGS  <- TRUE;
  
    if (is.null(genoud.parms$unif.seed))
      genoud.parms$unif.seed  <- 812821;
    
    if (is.null(genoud.parms$int.seed))
      genoud.parms$int.seed  <- 53058;
    
    if (is.null(genoud.parms$print.level))
      genoud.parms$print.level  <- 0;
    
    if (is.null(genoud.parms$share.type))
      genoud.parms$share.type  <- 0;
    
    if (is.null(genoud.parms$instance.number))
      genoud.parms$instance.number  <- 0;
    
    if (is.null(genoud.parms$output.path))
      genoud.parms$output.path  <- "stdout";
    
    if (is.null(genoud.parms$output.append))
      genoud.parms$output.append  <- FALSE;
    
    if (is.null(genoud.parms$project.path))
      genoud.parms$project.path  <- "/dev/null";
    
    if (is.null(genoud.parms$P1))
      genoud.parms$P1  <- 50;
    
    if (is.null(genoud.parms$P2))
      genoud.parms$P2  <- 50;
    
    if (is.null(genoud.parms$P3))
      genoud.parms$P3  <- 50;
    
    if (is.null(genoud.parms$P4))
      genoud.parms$P4  <- 50;
    
    if (is.null(genoud.parms$P5))
      genoud.parms$P5  <- 50;
    
    if (is.null(genoud.parms$P6))
      genoud.parms$P6  <- 50;
    
    if (is.null(genoud.parms$P7))
      genoud.parms$P7  <- 50;
    
    if (is.null(genoud.parms$P8))
      genoud.parms$P8  <- 50;
    
    if (is.null(genoud.parms$P9))
      genoud.parms$P9  <- 0  ;

    return(genoud.parms);
  } #end genoudParms


genoudRob <- function(fn,nvars,starting.values,genoud.parms)
{
  #set static defaults
  max  <- FALSE
  gradient.check  <- FALSE
  data.type.int  <- FALSE
  hessian  <- FALSE  
  roptim <- TRUE;

  #load up genoud.parms
  pop.size  <- genoud.parms$pop.size;
  max.generations  <- genoud.parms$max.generations;
  wait.generations  <- genoud.parms$wait.generations;
  hard.generation.limit  <- genoud.parms$hard.generation.limit;
  MemoryMatrix  <- genoud.parms$MemoryMatrix;
  Debug  <- genoud.parms$Debug;
  Domains  <- genoud.parms$Domains;
  scale.domains  <- genoud.parms$scale.domains;
  boundary.enforcement  <- genoud.parms$boundary.enforcement;
  solution.tolerance  <- genoud.parms$solution.tolerance;
  BFGS  <- genoud.parms$BFGS;
  unif.seed  <- genoud.parms$unif.seed;
  int.seed  <- genoud.parms$int.seed;
  print.level  <- genoud.parms$print.level;
  share.type  <- genoud.parms$share.type;
  instance.number  <- genoud.parms$instance.number;
  output.path  <- genoud.parms$output.path;
  output.append  <- genoud.parms$output.append;
  project.path  <- genoud.parms$project.path;
  P1  <- genoud.parms$P1;
  P2  <- genoud.parms$P2;
  P3  <- genoud.parms$P3;
  P4  <- genoud.parms$P4;
  P5  <- genoud.parms$P5;
  P6  <- genoud.parms$P6;
  P7  <- genoud.parms$P7;
  P8  <- genoud.parms$P8;
  P9  <- genoud.parms$P9;

  #we always have starting, but leave this check in. 
  #do we have starting values?
  if (is.null(starting.values)) {
    nStartingValues <- 0;
    parm.vec  <- rep(1,nvars)
  }
  else {
    nStartingValues <- 1;
    parm.vec  <- starting.values;
  }

  # let's create the Domains if none have been passed.
  if (!(is.matrix(Domains)))
    {
      Domains <- matrix(nrow=nvars, ncol=2);
      for (i in 1:nvars)
        {
          Domains[i,1] <- parm.vec[i] - abs(parm.vec[i])*scale.domains;
          Domains[i,2] <- parm.vec[i] + abs(parm.vec[i])*scale.domains;
        } # end of for loop
    } # end of Domains if
  
  #MemoryMatrix
  if (is.null(MemoryMatrix)) {
    MemoryMatrix <- TRUE;

    if (nvars > 20) {
      if (print.level > 0) {
        cat("\nWARNING: Since the number of parameters is greater than 20,\nWARNING: MemoryMatrix has been turned off by default.\nWARNING: You may turn it on using the MemoryMatrix flag.\nWARNING: This option increases speed at the cost of extra memory usage.\n\n")
        MemoryMatrix <- FALSE;
      }
    }
  }

  #set output.type
  if (output.path=="stdout")
    {
      output.type <- 0;
    }
  else
    {
      if (output.append)
        {
          output.type <- 2;
        }
      else
        {
          output.type <- 1;
        }
    }

  # create the P vector
  P <- vector(length=9, mode="numeric");
  P[1] <- P1; P[2] <- P2; P[3] <- P3; P[4] <- P4;
  P[5] <- P5; P[6] <- P6; P[7] <- P7; P[8] <- P8;
  P[9] <- P9;

  # has the user provided any seeds?
  if (unif.seed==812821 && int.seed==53058)
    provide.seeds <- FALSE
  else
    provide.seeds <- TRUE;

  if (max==FALSE)
        {
          g.scale <- 1;
        }
  else
    {
      g.scale <- -1;
    }

  #optim st
  genoud.optim.wrapper101 <- function(foo.vals)
    {
      ret <- optim(foo.vals, fn=as.function(fn), method="BFGS",
                  control=list(fnscale=g.scale));
      return(c(ret$value,ret$par));
    } # end of genoud.optim.wrapper101


  gout <- .Call("rgenoud", as.function(fn), new.env(),
               as.integer(nvars), as.integer(pop.size), as.integer(max.generations),
               as.integer(wait.generations),
               as.integer(nStartingValues), as.vector(starting.values),
               as.vector(P), as.matrix(Domains),
               as.integer(max), as.integer(gradient.check), as.integer(boundary.enforcement),
               as.double(solution.tolerance), as.integer(BFGS), as.integer(data.type.int),
               as.integer(provide.seeds), as.integer(unif.seed), as.integer(int.seed),
               as.integer(print.level), as.integer(share.type), as.integer(instance.number),
               as.integer(MemoryMatrix), as.integer(Debug),
               as.character(output.path), as.integer(output.type), as.character(project.path),
               as.integer(hard.generation.limit),
               as.function(genoud.optim.wrapper101), as.integer(roptim));

  if (hessian==TRUE)
    {
      hes <- optim(gout[5:(nvars+4)], fn, method="BFGS", hessian=TRUE,
                  control=list(fnscale=g.scale));
      
      hes <- hes$hessian;

      ret <- list(value=gout[1], generations=gout[2], peakgeneration=gout[3], popsize=gout[4],
                 par=gout[5:(nvars+4)], gradients=gout[(nvars+5):(nvars+nvars+4)],
                 operators=gout[(nvars+nvars+5):(nvars+nvars+9+4)],
                 hessian=hes);
    }
  else
    {
      ret <- list(value=gout[1], generations=gout[2], peakgeneration=gout[3], popsize=gout[4],
                 par=gout[5:(nvars+4)], gradients=gout[(nvars+5):(nvars+nvars+4)],
                 operators=gout[(nvars+nvars+5):(nvars+nvars+9+4)]);
    }

  return(ret);
} #end of genoudRob()


    
