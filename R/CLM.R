TraceEfficientSets<-function(C, mu, A, b, constraint_direction) {
  library(lpSolve)
  library(hash)
  corner_portfolios<-list()
  
  cat("Given input:","\n")
  cat("\nC","\n")
  print(C)
  cat("\nmu","\n")
  print(mu)
  cat("\nA","\n")
  print(A)
  cat("\nb","\n")
  print(b)
  cat("\nconstraint_direction","\n")
  print(constraint_direction)
  
  nassets<-dim(mu)[1]
  cat("Number of assets=",nassets,"\n")
  powers_of_two <- 2^(0:(nassets - 1)) #Idea: https://spartanideas.msu.edu/2019/02/21/binary-integer-conversion-in-r/
  nspaces<-(2^nassets-1)
  cat("Number of (portfolio) spaces=",nspaces,"\n")
  travelled_spaces<-hash()
  cat("\nRequired variables:","\n")
  M<-rbind(cbind(C, t(A) ), cbind(A, matrix(0, dim(A)[1], dim(A)[1])))
  cat("\nM","\n")
  print(M)
  R<-rbind(matrix(0, dim(A)[2]), b)
  cat("\nR","\n")
  print(R)
  S<-rbind(cbind(mu), matrix(0, dim(A)[1]))
  dimnames(S)<-NULL
  cat("\nS","\n")
  print(S)
  
  cat("\nFinding the portfolio with maximum E","\n")
  X1<-lp(objective.in = mu,
         const.mat = A,
         const.rhs = b,
         const.dir = constraint_direction,
         direction = "max")
  print(X1)
  cps<-ifelse(X1$solution==0,0,1)
  travelled_spaces<-SetTravelled(travelled_spaces,cps, powers_of_two)
  cat("Portfolio with maximum E (",round(t(as.matrix(X1$solution))%*%mu,3),") is (",round(X1$solution,2),") and it resides in portfolio space (",cps,")","\n\n")
  corner_portfolios[[length(corner_portfolios)+1]]<-cps
  cat("\nTracing efficient segements along critical lines residing in possible portfolio spaces","\n\n")
  round<-1
  repeat {
    N<-SetUnitCrosses(M, cps)
    N<-solve(N)
    N<-SetZeroCrosses(N, cps)
    
    T<-N%*%R
    U<-N%*%S
    cat("\nCritical line",round," in portfolio space (",cps,")","\n")
    cat("\nT","\n")
    print(T)
    cat("\nU","\n\n")
    print(U)
    
    if (exists("lambda_E")) {
      previous_lambda_E<-lambda_E
    }
    
    cat("\nCritical line",round,"intersections with other critical lines:","\n")
    ple<-0
    if (exists("lambda_E")) {
      ple<-lambda_E
    }
    variables<-GetNextSegment(mu,M,U,T,cps,travelled_spaces,ple, powers_of_two)
    add_variable<-variables$add_variable
    remove_variable<-variables$remove_variable
    le<-variables$lambda_E
    if (le<=0) {
      cat("\nThe points in portfolio space (",cps,") are efficient for lambda_E from",previous_lambda_E,"to 0","\n")
      lep<-NULL
      for (lambda in seq(round(as.numeric(previous_lambda_E),2), 0, -0.01)) {
        efficient_set<-round(head(T + lambda*U,nassets),2)
        if (any(efficient_set<0)) {
          break
        } else {
          lep<-efficient_set
        }
      }
      if (is.null(lep) == FALSE) {
        cat("Portfolio with minimum V (",t(lep)%*%C%*%lep,") is (",lep,"). E is",t(lep)%*%mu,"\n")
        corner_portfolios[[length(corner_portfolios)+1]]<-lep
        cat("Problem is over","\n")
        cat("A\n")
        break  
      }
      else {
        lep<-c(rep(0,nassets-1),1)
        cat("Portfolio with minimum V (",t(lep)%*%C%*%lep,") is (",lep,"). E is",t(lep)%*%mu,"\n")
        corner_portfolios[[length(corner_portfolios)+1]]<-lep
        cat("Problem is over","\n")
        cat("B\n")
        break
      }
    }
    lambda_E<-variables$lambda_E
    
    corner_portfolio<-NULL
    corner_average<-0
    if (exists("previous_lambda_E")) {
      cat("\nThe points in portfolio space (",cps,") are efficient for lambda_E from",previous_lambda_E,"to ",lambda_E,"where critical line",round,"intercepts critical line",round+1,"\n")
      for (lambda in seq(round(as.numeric(previous_lambda_E),2), round(as.numeric(lambda_E),2), -0.02)) {
        efficient_set<-round(head(T + lambda*U,nassets),2)
      }
      corner_portfolio<-round(head(T+U%*%lambda_E,nassets),2)
      corner_average<-t(head(T+U%*%lambda_E,nassets))%*%mu
    }
    
    if (add_variable>0) {
      pps<-cps
      cps[add_variable]<-1
      cat("\nNext IN variable is",add_variable,"and so next efficient set segment lies along critical line",round+1,"in portfolio space (", cps,")","\n")
      if (is.null(corner_portfolio)==FALSE) {
        cat("\nThe corner portfolio at which the efficient set turns from portfolio space (",pps,") to (",cps,") is (",corner_portfolio,") with E",round(corner_average,3),"\n")
        corner_portfolios[[length(corner_portfolios)+1]]<-corner_portfolio
      }
      travelled_spaces<-SetTravelled(travelled_spaces,cps, powers_of_two)
    } else if (remove_variable>0) {
      pps<-cps
      cps[remove_variable]<-0
      cat("\nNext OUT variable is",remove_variable,"and so next efficient set segment lies along critical line",round+1,"in portfolio space (", cps,")","\n")
      if (is.null(corner_portfolio)==FALSE) {
        cat("\nThe corner portfolio at which the efficient set turns from portfolio space (",pps,") to (",cps,") is (",corner_portfolio,") with E",round(corner_average,3),"\n")  
        corner_portfolios[[length(corner_portfolios)+1]]<-corner_portfolio
      }
      travelled_spaces<-SetTravelled(travelled_spaces,cps, powers_of_two)
    }
    round<-round+1
  }
  return(corner_portfolios)
}

SetZeroCrosses<-function(M, inout) {
  zeroRowCol<-rep(0,dim(M)[1]) #M is symmetric
  psi<-1
  for (isOut in inout==0) {
    if (isOut==TRUE) {
      M[,psi]<-zeroRowCol
      M[psi,]<-zeroRowCol
    }
    psi=psi+1
  }
  return(M)
}

SetUnitCrosses<-function(M, inout) {
  M<-SetZeroCrosses(M,inout)
  psi<-1
  for (isOut in inout==0) {
    if (isOut==TRUE) {
      M[psi,psi]=1
    }
    psi=psi+1
  }
  return(M)
}

IsTravelled<-function(travelled_spaces, ps, powers_of_two) {
  key<-as.character(ps %*% powers_of_two)
  if (is.null(travelled_spaces[[key]]) || travelled_spaces[[key]] == FALSE) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

SetTravelled<-function(travelled_spaces, ps, powers_of_two) {
  key<-as.character(ps %*% powers_of_two)
  travelled_spaces[[key]]<-TRUE
  return(travelled_spaces)
}

GetNextSegment<-function(mu, M, U, T, cps, travelled_spaces, previous_lambda_E, powers_of_two) {
  add_variable<-0
  remove_variable<-0
  largest_lambda_E<-0
  psi<-1
  for (isOut in cps==0) {
    if (isOut==TRUE) {
      ps<-cps
      ps[psi]<-1
      travelled<-IsTravelled(travelled_spaces, ps, powers_of_two)
      if (travelled==TRUE) {
        cat("If added variable",psi,"then next critical line would be in portfolio space (",ps,"), but found that it is already travelled","\n")
      } else {
        l<-solve(mu[psi]-M[,psi]%*%U, M[,psi]%*%T)
        cat("If added variable",psi,"then next critical line would be in portfolio space (",ps,") and lambda_E is",l,"\n")
        if (l>largest_lambda_E) {
          if (previous_lambda_E>0 && l>previous_lambda_E) {
            cat("   ....but that is greater than the previous one","\n")
          } else {
            largest_lambda_E<-l
            add_variable<-psi
            remove_variable<-0
          }
        }
      }
    } else {
      if (sum(cps)>1) {
        ps<-cps
        ps[psi]<-0
        travelled<-IsTravelled(travelled_spaces, ps, powers_of_two)
        if (travelled==TRUE) {
          cat("If removed variable",psi,"then next critical line would be in portfolio space (",ps,") but found that it is already travelled","\n")
        } else {
          l<-solve(U[psi],-T[psi])
          cat("If removed variable",psi,"then next critical line would be in portfolio space (",ps,") and lambda_E is",l,"\n")
          if (l>largest_lambda_E) {
            if (previous_lambda_E>0 && l>previous_lambda_E) {
              cat("   ....but that is greater than the previous one","\n")
            } else {
              largest_lambda_E<-l
              remove_variable<-psi
              add_variable<-0
            }
          }
        }
      }
    }
    psi=psi+1
  }
  return(list(lambda_E=largest_lambda_E, add_variable=add_variable, remove_variable=remove_variable))
}
