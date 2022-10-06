#!/usr/local/R/3.6.2/bin/Rscript
####Numerical experiment

####TV constarned Max-flow  
####with JD data

# Load linear programming library
# install.packages("lpSolve")
library("lpSolve")

params = read.csv("params.csv", header = TRUE)
n1 <- params$n1
n2 <- params$n2
n <- params$n
l <- params$l
k <- params$k
print(c(n1,n2,n,k,l))
print(l*n+k^2)
p <- read.csv("priceset_std.csv", header = TRUE)
p <- unlist(p[,1])

l=length(p)
print(l)
## Index of the dataset
NO = 4
# # generate feature x
features1 <- c('user_level', 'age', 'marital_status', 'education', 'city_level', 'purchase_power', 'gender')
X = read.csv("X_std_train.csv", header = TRUE)
print(str(X))
X = X[, features1]
print(str(X))
X = as.matrix(X)
print(dim(X))
print(str(X))

#### for delta in drange
# drange <- c()
# for (k in 1:4) {
#   i <- 5-k
#   drange[i] <- 10^(-k-1)
# }
# drange <- c(c(10^(-6), 10^(-3)), seq(from = 0.01, to = 0.1, length.out = 7), seq(from = 0.15, to = 1, length.out = 11))
# length(drange) # 20
drange <- seq(from = 0.05, to = 1, by = 0.05)

# column_values = ['intercept', 'final_unit_price', 'user_level', 'age', 'marital_status', 'education', 'city_level', 'purchase_power', 'gender']
features2 <- c('intercept', 'user_level', 'age', 'marital_status', 'education', 'city_level', 'purchase_power', 'gender')
a_df <- read.csv("reg.coef.csv", header = TRUE)
a <- a_df[, features2]
a <- as.matrix(a)
print(a)
print(a_df[,"final_unit_price"])
            # Simulate V and compute Surplus
            print(dim(t(X)))
           # print(dim(matrix(rbind(rep(1,n),t(X)), nrow = 1+length(features1))))
            ax <- matrix(a, nrow = 1) %*% matrix(rbind(rep(1,n),t(X)), nrow = 1+length(features1))
            cx <- - ax/a_df[,"final_unit_price"]
            summary(t(cx)) # 2.002238 6.997109
            summary(t(ax)) # 3.004476 8.994217
            # Customer Valuation
            V <- c()
            for (i in 1:n) {
              # keep track of random seed
              set.seed(i)
              # Linear demand - Uniform distribution
              V[i] <- runif(1,0,cx[i]) 
            }
            summary(V)
            # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
            # 0.000236 0.950737 1.918675 2.203933 3.148514 6.826074
            # Demand
    # create coefficient vector
    d <- read.csv("demand_allp.csv", header = TRUE)
    d_vec <- read.csv("demand_vec_allp.csv", header = TRUE)
    d_vec <- unlist(d_vec[,1])
    print("d",d_vec)
            #Emperical Revenue 
  # r <- matrix(ncol = l, nrow = n)
  # obj_coef <- c()
r <- read.csv("revenue_allp.csv", header = TRUE)
r <- as.matrix(r)
 obj_coef <- read.csv("obj_coef.csv", header = TRUE)
 obj_coef <- unlist(obj_coef[,1])
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
  # initialize counter for delta
  m <- 0
  opt_res <- matrix(nrow = length(drange), ncol = 10+k*l+l+1+l+1)
  colnames(opt_res) <- c("delta", "opt_obj", "opt_obj1","opt_obj2","opt_surplus", "opt_surplus1","opt_surplus2","opt_welfare","opt_welfare1","opt_welfare2",paste0('alp1_',seq(1,l,1)),paste0('alp2_',seq(1,l,1)),paste0('alp12_diff_',seq(1,l,1)),'alp_diff',paste0('gini_coef',seq(1,l,1)),'gini')
  for (delta in drange) {
    m <- m+1
    f.obj <- c(obj_coef, rep(0,k*l+k*(k-1)*l/2))
    f.con <- mat
    f.dir <- c(rep("<=",k*l), rep("<=", n), rep("<=", k*l), "<=", rep("=",k), rep(">=", l*n + k*l+k*(k-1)*l/2))
    f.rhs <- c(rep(0,k*l), rep(1, n), rep(0, k*l), delta, rep(1,k), rep(0, l*n + k*l+k*(k-1)*l/2))
    
    # Results
    opt_obj <- lp ("max", f.obj, f.con, f.dir, f.rhs)$objval
    opt_sol <- lp ("max", f.obj, f.con, f.dir, f.rhs)$solution
    print(c("optimal obj",m, ":", opt_obj))
    # Compute the optimal price distribution difference
    opt_sol_alp_diff <- abs(opt_sol[(l*n+1):(l*n+l)] - opt_sol[(l*n+l+1):(l*n+k*l)])
    opt_alp_diff <- sum(opt_sol_alp_diff)
    
    # Compute Gini coefficient for each price and Gini index
    gini_coef <- abs(opt_sol[(l*n+1):(l*n+l)] - opt_sol[(l*n+l+1):(l*n+k*l)]) / (k*(opt_sol[(l*n+1):(l*n+l)] + opt_sol[(l*n+l+1):(l*n+k*l)]))
    gini <- sum(gini_coef, na.rm = TRUE)
     
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
    opt_surplus <- opt_sol[1:(l*n)] * d_vec * pmax(V - rep(p, each = n),0)
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
      opt_res[m,] <- c(delta, opt_obj,opt_obj1,opt_obj2, opt_surplus,opt_surplus1,opt_surplus2, opt_welfare,opt_welfare1,opt_welfare2,opt_sol[(l*n+1):(l*n+k*l)],opt_sol_alp_diff,opt_alp_diff,gini_coef,gini)
    }

write.csv(opt_res, file = paste0("res_TVconst_VD_multi-delta(",min(drange), max(drange),")_JDdata", NO,".csv"))

#####################sanity Check
warnings()
