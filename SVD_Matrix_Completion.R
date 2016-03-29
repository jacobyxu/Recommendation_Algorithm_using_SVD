### The function to get raw matrix M
Get_M_NA <- function (u){
  M <- matrix(NA, max(u[,1]),max(u[,2]))
  for (i in 1:nrow(u)){
    x = u[i,1]
    y = u[i,2]
    M[x,y] = u[i,3]
  }
  return(M)
}
### The function of SoftImpute Algorithm by SVD
Self_SVD <- function (M, lambda = 0, maxiteN = 100, epsilon = 1e-5 ){
  library(svd)
  # record the all plot of reliable values
  N <- matrix(NA,2,nrow(M)*ncol(M))
  t <- 1
  for (i in 1: nrow(M)){
    for (j in 1:ncol(M)){
      if (!is.na(M[i,j])){
        N[1,t] <- i
        N[2,t] <- j
        t <- t+1
      }
    }
  }
  N <- N[,1:t]
  Itr <- c(N[1,])
  Jtr <- c(N[2,])
  # Change NA to 0
  for (i in 1: nrow(M)){
    for (j in 1:ncol(M)){
      if (is.na(M[i,j]))  M[i,j] <- 0
    }
  }
  # Algorithm
  ## set default Z_old, Z_new
  Z_old <- matrix(0, nrow(M),ncol(M))
  Z_new <- matrix(0, nrow(M),ncol(M))
  ## Repeat
  count_ite <- 0
  for (j in 1:maxiteN) {
    ### Caculate Z_new
    for (i in 1:length(Itr)) Z_old[Itr[i],Jtr[i]] <- 0
    Z_temp <- Z_old + M
    SVD_Z_temp <- svd(Z_temp)
    SVD_Z_temp$d <- max(0, SVD_Z_temp$d - lambda)
    d <- matrix(0, ncol(SVD_Z_temp$u), ncol(SVD_Z_temp$v))
    for (k in 1:length(SVD_Z_temp$d)) d[k,k] <- SVD_Z_temp$d[k]
    Z_new <- SVD_Z_temp$u %*% d %*% t(SVD_Z_temp$v)
    ### check the eps
    current_eps <- (norm(Z_new - Z_old, type = "F"))/(norm(Z_old, type = "F")+0.001)
    if (current_eps < epsilon) break
    ### Assign
    Z_old <- Z_new
    count_ite <- count_ite + 1 # count iterations
  }
  ## Assign
  for (i in 1:length(Itr)) Z_old[Itr[i],Jtr[i]] <- 0
  Z <- Z_old + M
  out = list(Z,count_ite, current_eps)
  return(out)
}
# Get the raw matrix M
t0 <- Sys.time() # start the timer
u <- read.table("u.data")
u <- u[,1:3]
M <- Get_M_NA(u)
Result <- Self_SVD(M, lambda = 5, maxiteN = 500, epsilon = 0.1)
Z <- Result[[1]]
iteTimes <- Result[[2]]
eps <- Result[[3]]
write.csv(Z,"Z.csv")
t1 <- Sys.time() # stop the timer
running_time <- t1 - t0
iteTimes
running_time
eps





# Seperate to training and testing sets
training_set_nrow <- sample(nrow(M),floor(nrow(M)/(sqrt(3)+1)*sqrt(3)))
training_set_ncol <- sample(ncol(M),floor(ncol(M)/(sqrt(3)+1)*sqrt(3)))
training_set <- M[training_set_nrow,]
training_set <- training_set[,training_set_ncol]
N_row <- matrix(1:nrow(M),1,nrow(M))
N_col <- matrix(1:ncol(M),1,ncol(M))
testing_set_nrow <- c(N_row[1,-training_set_nrow])
testing_set_ncol <- c(N_col[1,-training_set_ncol])
testing_set <- M[testing_set_nrow,]
testing_set <- testing_set[,testing_set_ncol]



write.csv(training_set,"Training.csv")
write.csv(testing_set,"Testing.csv")
t1 <- Sys.time() # stop the timer
running_time <- t1 - t0












