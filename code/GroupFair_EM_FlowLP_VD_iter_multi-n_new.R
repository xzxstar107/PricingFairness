#############
#!/usr/local/R/3.6.2/bin/Rscript
####Multi-n experiment

#############
####EM constarned Max-flow  
####with simulated data
# Load linear programming library
library(lpSolve)

# parameter intialization
seed0 <- 2020
set.seed(seed0)

# parameter intialization
# delta <- 10^(-5)
delta <- 10^(-1)
l <- 10
k <- 2
pmin <- 1
pmax <- 3
p <- seq(length.out =l,from = pmin, to = pmax) 
# n_param0 <- 2
# paramvec0 <- seq(length.out =n_param0,from = 1, to = 2)
# n_param <- 3
# paramvec <- seq(length.out =n_param,from = 1, to = 5)
paramind <- 1
c <- c(1, 1, 1)
a <- c(1, 1, 1)

#### for n in nrange
n <- 2000
nrange <- seq(50, n/2, by = 50)
##
opt_res <- matrix(ncol = k+3*(k+1)+k*l , nrow = length(nrange))
colnames(opt_res) <- c("n1","n2", "opt_obj", "opt_obj1","opt_obj2", "opt_surplus","opt_surplus1","opt_surplus2", 
                       "opt_welfare","opt_welfare1","opt_welfare2", rep("alp1",l), rep("alp2",l))
params_init <- c(a, c, delta, l, k, n, p)

# Multi-n loop 
ind <- 0
for (n1 in nrange) {
  ind <- ind+1
  n2 <- n-n1
  set.seed(2020)
  # generate feature x
  x <- runif(n1, 0, 1)
  s <- c()
  s[1:n1] <- 0
  # generate feature y
  y <- runif(n2, 0, 1)
  s[(n1+1): n] <- 1
  
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
  #Emperical Revenue 
  r <- matrix(ncol = l, nrow = n)
  obj_coef <- c()
  ##
  for(j in 1:l){
    r[, j] <- ax * (1-p[j]/cx) * p[j]
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
  mat <- matrix(0, ncol = l*n + k*l+l^2, nrow = k*l+n+1+k*l+k+n*l+k*l+l^2)
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
  # EM distance
  # distane matrix
  dmat <- matrix(nrow = l, ncol = l)
  for (i in 1:l) {
    for (j in 1:l){
      dmat[i,j] <- abs(p[i] - p[j])
      dmat[j,i] <- dmat[i,j]
    }
  }
  # coefficient: distance vector 
  dvec <- c()
  for (i in 1:l) {
    dvec <- c(dvec, dmat[i,])
  }
  #Normalized 
  # Fairness constraint
  st_ind <- k*l+n
  mat[st_ind+1, (l*n+k*l+1):(l*n+k*l+l^2)] <- dvec/(pmax - pmin)
  # \alpha1
  st_ind <- k*l+n+1
  for (i in 1:l) {
    mat[st_ind+i, l*n+l+i] <- -1
  }
  # hij
  for (i in 1:l) {
    for (j in 1:l) {
      mat[st_ind+i, (l*n+k*l+i+(j-1)*l)] <- 1
    }
  }
  #\alpha2
  st_ind <- k*l+n+1+l
  for (i in 1:l) {
    mat[st_ind+i, l*n+i] <- -1
  }
  #hij
  for (i in 1:l) {
    mat[st_ind+i, (l*n+k*l+1+(i-1)*l):(l*n+k*l+l+(i-1)*l)] <- rep(1,l)
  }
  #\alphas
  st_ind <- k*l+n+1+k*l
  for (i in 1:k) {
    mat[st_ind+i, (l*n+(i-1)*l+1):(l*n+(i-1)*l+l)] <- rep(1, l)
  }
  
  # Variables nonnegative condition
  st_ind <- k*l+n+1+k*l+k
  mat[(st_ind+1):(st_ind+l*n + k*l+l^2), ] <- diag(x=1, nrow = l*n + k*l+l^2, ncol = l*n + k*l+l^2)
  # head(mat)
  # dim(mat) 
  
  f.obj <- c(obj_coef, rep(0, k*l+l^2))
  f.con <- mat
  f.dir <- c(rep("<=",k*l), rep("<=", n), "<=", rep("=", k*l), rep("=",k), rep(">=", l*n + k*l+l^2))
  f.rhs <- c(rep(0,k*l), rep(1, n), delta, rep(0, k*l), rep(1,k), rep(0, l*n + k*l+l^2))
  
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
  opt_surplus <- opt_sol[1:(l*n)] * max(V - rep(p, each = n),0)
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
  opt_res[ind,] <- c(n1, n2, opt_obj, opt_obj1,opt_obj2, opt_surplus,opt_surplus1,opt_surplus2, 
                     opt_welfare,opt_welfare1,opt_welfare2, opt_sol[(l*n+1): (l*n+k*l)])
}

# Output settings and results
write.csv(params_init, file = paste0('paramsSim_EM_n_new.csv'))
write.csv(opt_res, file = paste0('res_EM_n_',paramind,'_new.csv'))



