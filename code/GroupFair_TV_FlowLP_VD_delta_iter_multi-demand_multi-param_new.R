#!/usr/local/R/3.6.2/bin/Rscript
####Multi-parameter experiment

####TV constarned Max-flow 
####under multiple demand models 
####with simulated data

# Load linear programming library
# install.packages("lpSolve")
library("lpSolve")

# initial set-up
n1 <- 1000
n2 <- 1000
n <- n1 + n2
seed0 <- 2020
set.seed(seed0)
# generate feature x
x <- runif(n, 0, 1)
s <- c()
s[1:n1] <- 0
s[(n1+1): n] <- 1
# length(s)
# summary(x)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0007724 0.2442752 0.4872140 0.4958653 0.7439779 0.9999953  
# summary(s)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0     0.0     0.5     0.5     1.0     1.0 
# parameter intialization
# delta <- 10^(-5)
# delta <- 10^(-1)
l <- 10
k <- 2
pmin <- 1
pmax <- 3
p <- seq(length.out =l,from = pmin, to = pmax) 
n_param0 <- 2
paramvec0 <- seq(length.out =n_param0,from = 1, to = 2)
n_param <- 2
paramvec <- seq(length.out =n_param,from = 1, to = 2)

#### for delta in drange
# drange <- c()
# for (k in 1:4) {
#   i <- 5-k
#   drange[i] <- 10^(-k-1)
# }
# drange <- c(c(10^(-6), 10^(-3)), seq(from = 0.01, to = 0.1, length.out = 7), seq(from = 0.15, to = 1, length.out = 11))
# length(drange) # 20
drange <- seq(from = 0.05, to = 1, by = 0.05)

paramind <- 0
# Fix constant
i1 <- 1
j1 <- 1
## Iteration over multi-parameters
for (i2 in paramvec0) {
  for (i3 in paramvec0) {
      for (j2 in paramvec) {
        for (j3 in paramvec) {
            # sample parameters
            paramind <- paramind + 1
            c <- c(i1, i2, i3)
            a <- c(j1, j2, j3)
            # Simulate V and compute Surplus
            ax <- matrix(a, nrow = 1) %*% matrix(rbind(rep(1,n),x,s), nrow = 3)
            cx <- matrix(c, nrow = 1) %*% matrix(rbind(rep(1,n),x,s), nrow = 3)
            # summary(t(cx)) # 2.002238 6.997109
            # summary(t(ax)) # 3.004476 8.994217
            # Customer Valuation
            V <- c()
            for (i in 1:n) {
              # keep track of random seed
              set.seed(i)
              # Linear demand - Uniform distribution
              V[i] <- runif(1,0,cx[i]) 
            }
            # summary(V)
            # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
            # 0.000236 0.950737 1.918675 2.203933 3.148514 6.826074
            #Emperical Revenue with exponential demand 
  r <- matrix(ncol = l, nrow = n)
  obj_coef <- c()
  for(j in 1:l){
     r[, j] <- ax * exp(-cx*p[j]) * p[j]
    # create coefficient vector
    obj_coef <- c(obj_coef, r[, j])
  }
  # create coefficient vector for group 1
  obj_coef1 <- c() #l*n1
  for (j in 1:l) {
    obj_coef1 <- c(obj_coef1, obj_coef[((j-1)*n+1):((j-1)*n+n1)])
  }
  # create coefficient vector for group 2
  obj_coef2 <- c() #l*n2
  for (j in 1:l) {
    obj_coef2 <- c(obj_coef2, obj_coef[((j-1)*n+n1+1):((j-1)*n+n1+n2)])
  }
  
  # create coeffecients for constrants
  mat <- matrix(0, ncol = l*n + k*l+k*(k-1)*l/2, nrow = k*l+n+k*l+1+k+n*l+k*l+k*(k-1)*l/2)
  for (j in 1:l) {
    mat[j, ((j-1)*n+1):((j-1)*n+n1)] <- rep(1, n1)
    mat[l+j, ((j-1)*n+n1+1):((j-1)*n+n1+n2)] <- rep(1, n2)
    mat[j, l*n+j] <- -n1
    mat[l+j, l*n+l+j] <- -n2
  }
  
  st_ind <- k*l
  for (i in 1:n1) {
    for (j in 1:l) {
      mat[st_ind+i,(j-1)*n+i] <- 1
    }
  }
  
  st_ind <- k*l+n1
  for (i in 1:n2) {
    for (j in 1:l) {
      mat[st_ind+i,(j-1)*n+n1+i] <- 1
    }
  }
  
  # introduce instrumental variables for Fairness constraint
  st_ind <- k*l+n
  for (j in 1:l) {
    mat[st_ind+2*(j-1)+1, c(n*l+j, n*l+(k-1)*l+j, n*l+k*l+j)] <- c(1, -1, -1)
    mat[st_ind+2*(j-1)+2, c(n*l+j, n*l+(k-1)*l+j, n*l+k*l+j)] <- c(-1, 1, -1)
  }
  # Fairness constraint
  st_ind <- k*l+n+k*l
  mat[st_ind+1, (l*n+k*l+1):(l*n+k*l+l)] <- rep(1, l)
  # 
  st_ind <- k*l+n+k*l+1
  for (i in 1:k) {
    mat[st_ind+i, (l*n+(i-1)*l+1):(l*n+(i-1)*l+l)] <- rep(1, l)
  }
  
  # Variables nonnegative condition
  st_ind <- k*l+n+k*l+1+k
  mat[(st_ind+1):(st_ind+l*n + k*l+k*(k-1)*l/2), ] <- diag(x=1, nrow = l*n + k*l + k*(k-1)*l/2, 
                                                           ncol = l*n + k*l + k*(k-1)*l/2)
  # head(mat)
  # dim(mat) 
  
  ## iteration for delta
     opt_res <- matrix(nrow = length(drange), ncol = 6+10+k*l)
     alp_name <- c()
     for (i in 1:k) {
       alp_name <- c(alp_name, paste0("alp",i,"_",1:l))
     }
     colnames(opt_res) <- c("a1","a2","a3","c1","c2","c3","delta", "opt_obj", "opt_obj1","opt_obj2","opt_surplus", "opt_surplus1","opt_surplus2","opt_welfare","opt_welfare1","opt_welfare2", alp_name)
  # initialize counter for delta 
  m <- 0    
  for (delta in drange) {
    m <- m+1
    f.obj <- c(obj_coef, rep(0,k*l+k*(k-1)*l/2))
    f.con <- mat
    f.dir <- c(rep("<=",k*l), rep("<=", n), rep("<=", k*l), "<=", rep("=",k), rep(">=", l*n + k*l+k*(k-1)*l/2))
    f.rhs <- c(rep(0,k*l), rep(1, n), rep(0, k*l), delta, rep(1,k), rep(0, l*n + k*l+k*(k-1)*l/2))
    
    # Results
    opt_obj <- lp ("max", f.obj, f.con, f.dir, f.rhs)$objval
    opt_sol <- lp ("max", f.obj, f.con, f.dir, f.rhs)$solution
    
    # Extract optimal solution for group1
    opt_sol1 <- c() #l*n1
    for (j in 1:l) {
      opt_sol1 <- c(opt_sol1, opt_sol[((j-1)*n+1):((j-1)*n+n1)])
    }
    # Extract optimal solution for group2
    opt_sol2 <- c() #l*n2
    for (j in 1:l) {
      opt_sol2 <- c(opt_sol2, opt_sol[((j-1)*n+n1+1):((j-1)*n+n1+n2)])
    }
    # Compute optimal revenue for group1 & 2
    opt_obj1 <- sum(opt_sol1 * obj_coef1)
    opt_obj2 <- sum(opt_sol2 * obj_coef2)
    
    # Compute surplus 
    # opt_surplus <- opt_sol[1:(l*n)] * max(V - rep(p, each = n),0)
    pvec <- c()
    for (j in 1:l){
      pvec <- c(pvec, rep(p[j],n))
    }
    opt_surplus <- opt_sol[1:(l*n)] * 1/rep(cx,l)* rep(ax,l)*exp(-rep(cx,l)*pvec)
    opt_surplus1 <-c()
    for (j in 1:l) {
      opt_surplus1 <- c(opt_surplus1, opt_surplus[((j-1)*n+1):((j-1)*n+n1)])
    }
    opt_surplus1 <- sum(opt_surplus1)
    opt_surplus2 <-c()
    for (j in 1:l) {
      opt_surplus2 <- c(opt_surplus2, opt_surplus[((j-1)*n+n1+1):((j-1)*n+n1+n2)])
    }
    opt_surplus2 <- sum(opt_surplus2)
    # Compute welfare
    opt_surplus <- sum(opt_surplus)
    opt_welfare <- opt_surplus + opt_obj
    opt_welfare1 <- opt_surplus1 + opt_obj1
    opt_welfare2 <- opt_surplus2 + opt_obj2
# Save results
    opt_res[m,] <- c(i1, i2, i3, j1, j2, j3, delta, opt_obj,opt_obj1,opt_obj2, opt_surplus,opt_surplus1,opt_surplus2, opt_welfare,opt_welfare1,opt_welfare2, opt_sol[(n*l+1):(n*l+k*l)])
    }
    # output results
    write.csv(opt_res, file = paste0("res_TVconstr_VD_multi-delta_expDemand_params_",i1,i2,i3,j1,j2,j3,"_new.csv"))
      }  
          }
          
        }
        
    }
      

# Logistic Demand
# Fix constant
i1 <- 1
j1 <- 1
t1 <- 1
## Iteration over multi-parameters
for (i2 in paramvec0) {
  for (i3 in paramvec0) {
    for (j2 in paramvec) {
      for (j3 in paramvec) {
        for (t2 in paramvec/5) {
          for (t3 in paramvec/5) {
            # sample parameters
            paramind <- paramind + 1
            c <- c(i1, i2, i3)
            a <- c(j1, j2, j3)
            b <- c(t1, t2, t3)
            # Simulate V and compute Surplus
            ax <- matrix(a, nrow = 1) %*% matrix(rbind(rep(1,n),x,s), nrow = 3)
            cx <- matrix(c, nrow = 1) %*% matrix(rbind(rep(1,n),x,s), nrow = 3)
            bx <- exp(matrix(b, nrow = 1) %*% matrix(rbind(rep(1,n),x,s), nrow = 3))
            # summary(t(cx)) # 2.002238 6.997109
            # summary(t(ax)) # 3.004476 8.994217
            # Customer Valuation
            V <- c()
            for (i in 1:n) {
              # keep track of random seed
              set.seed(i)
              # Linear demand - Uniform distribution
              V[i] <- runif(1,0,cx[i]) 
            }
            # summary(V)
            # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
            # 0.000236 0.950737 1.918675 2.203933 3.148514 6.826074
            #Emperical Revenue with logistic demand 
  r <- matrix(ncol = l, nrow = n)
  obj_coef <- c()
  for(j in 1:l){
     r[, j] <- bx/(bx+exp(-cx*p[j])) * p[j]
    # create coefficient vector
    obj_coef <- c(obj_coef, r[, j])
  }
  # create coefficient vector for group 1
  obj_coef1 <- c() #l*n1
  for (j in 1:l) {
    obj_coef1 <- c(obj_coef1, obj_coef[((j-1)*n+1):((j-1)*n+n1)])
  }
  # create coefficient vector for group 2
  obj_coef2 <- c() #l*n2
  for (j in 1:l) {
    obj_coef2 <- c(obj_coef2, obj_coef[((j-1)*n+n1+1):((j-1)*n+n1+n2)])
  }
  
  # create coeffecients for constrants
  mat <- matrix(0, ncol = l*n + k*l+k*(k-1)*l/2, nrow = k*l+n+k*l+1+k+n*l+k*l+k*(k-1)*l/2)
  for (j in 1:l) {
    mat[j, ((j-1)*n+1):((j-1)*n+n1)] <- rep(1, n1)
    mat[l+j, ((j-1)*n+n1+1):((j-1)*n+n1+n2)] <- rep(1, n2)
    mat[j, l*n+j] <- -n1
    mat[l+j, l*n+l+j] <- -n2
  }
  
  st_ind <- k*l
  for (i in 1:n1) {
    for (j in 1:l) {
      mat[st_ind+i,(j-1)*n+i] <- 1
    }
  }
  
  st_ind <- k*l+n1
  for (i in 1:n2) {
    for (j in 1:l) {
      mat[st_ind+i,(j-1)*n+n1+i] <- 1
    }
  }
  
  # introduce instrumental variables for Fairness constraint
  st_ind <- k*l+n
  for (j in 1:l) {
    mat[st_ind+2*(j-1)+1, c(n*l+j, n*l+(k-1)*l+j, n*l+k*l+j)] <- c(1, -1, -1)
    mat[st_ind+2*(j-1)+2, c(n*l+j, n*l+(k-1)*l+j, n*l+k*l+j)] <- c(-1, 1, -1)
  }
  # Fairness constraint
  st_ind <- k*l+n+k*l
  mat[st_ind+1, (l*n+k*l+1):(l*n+k*l+l)] <- rep(1, l)
  # 
  st_ind <- k*l+n+k*l+1
  for (i in 1:k) {
    mat[st_ind+i, (l*n+(i-1)*l+1):(l*n+(i-1)*l+l)] <- rep(1, l)
  }
  
  # Variables nonnegative condition
  st_ind <- k*l+n+k*l+1+k
  mat[(st_ind+1):(st_ind+l*n + k*l+k*(k-1)*l/2), ] <- diag(x=1, nrow = l*n + k*l + k*(k-1)*l/2, 
                                                           ncol = l*n + k*l + k*(k-1)*l/2)
  # head(mat)
  # dim(mat) 
  
  ## iteration for delta
     opt_res <- matrix(nrow = length(drange), ncol = 9+10+k*l)
     alp_name <- c()
     for (i in 1:k) {
       alp_name <- c(alp_name, paste0("alp",i,"_",1:l))
     }
     colnames(opt_res) <- c("a1","a2","a3","c1","c2","c3","b1","b2","b3","delta", "opt_obj", "opt_obj1","opt_obj2","opt_surplus", "opt_surplus1","opt_surplus2","opt_welfare","opt_welfare1","opt_welfare2", alp_name)
  # initialize counter for delta 
  m <- 0    
  for (delta in drange) {
    m <- m+1
    f.obj <- c(obj_coef, rep(0,k*l+k*(k-1)*l/2))
    f.con <- mat
    f.dir <- c(rep("<=",k*l), rep("<=", n), rep("<=", k*l), "<=", rep("=",k), rep(">=", l*n + k*l+k*(k-1)*l/2))
    f.rhs <- c(rep(0,k*l), rep(1, n), rep(0, k*l), delta, rep(1,k), rep(0, l*n + k*l+k*(k-1)*l/2))
    
    # Results
    opt_obj <- lp ("max", f.obj, f.con, f.dir, f.rhs)$objval
    opt_sol <- lp ("max", f.obj, f.con, f.dir, f.rhs)$solution
    
    # Extract optimal solution for group1
    opt_sol1 <- c() #l*n1
    for (j in 1:l) {
      opt_sol1 <- c(opt_sol1, opt_sol[((j-1)*n+1):((j-1)*n+n1)])
    }
    # Extract optimal solution for group2
    opt_sol2 <- c() #l*n2
    for (j in 1:l) {
      opt_sol2 <- c(opt_sol2, opt_sol[((j-1)*n+n1+1):((j-1)*n+n1+n2)])
    }
    # Compute optimal revenue for group1 & 2
    opt_obj1 <- sum(opt_sol1 * obj_coef1)
    opt_obj2 <- sum(opt_sol2 * obj_coef2)
    
    # Compute surplus 
    # opt_surplus <- opt_sol[1:(l*n)] * max(V - rep(p, each = n),0)
    pvec <- c()
    for (j in 1:l){
      pvec <- c(pvec, rep(p[j],n))
    }
    opt_surplus <- opt_sol[1:(l*n)] * 1/rep(cx,l)* log(rep(1,n*l)+rep(bx,l)*exp(-rep(cx,l)*pvec))
    opt_surplus1 <-c()
    for (j in 1:l) {
      opt_surplus1 <- c(opt_surplus1, opt_surplus[((j-1)*n+1):((j-1)*n+n1)])
    }
    opt_surplus1 <- sum(opt_surplus1)
    opt_surplus2 <-c()
    for (j in 1:l) {
      opt_surplus2 <- c(opt_surplus2, opt_surplus[((j-1)*n+n1+1):((j-1)*n+n1+n2)])
    }
    opt_surplus2 <- sum(opt_surplus2)
    # Compute welfare
    opt_surplus <- sum(opt_surplus)
    opt_welfare <- opt_surplus + opt_obj
    opt_welfare1 <- opt_surplus1 + opt_obj1
    opt_welfare2 <- opt_surplus2 + opt_obj2
# Save results
    opt_res[m,] <- c(i1, i2, i3, j1, j2, j3, t1,t2,t3, delta, opt_obj,opt_obj1,opt_obj2, opt_surplus,opt_surplus1,opt_surplus2, opt_welfare,opt_welfare1,opt_welfare2, opt_sol[(n*l+1):(n*l+k*l)])
    }
    # output results
    write.csv(opt_res, file = paste0("res_TVconstr_VD_multi-delta_logisticDemand_params_",i1,i2,i3,j1,j2,j3,t1,t2,t3"_new.csv"))
      }  
          }
          
        }
        
    }




