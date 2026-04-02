#' Heterogeneous Multiple-Inflated Negative Binomial (HMINB) Model
#'
#' Fits a multiple-inflated negative binomial model with covariate-dependent
#' inflation probabilities, using cumulative logit and fused LASSO penalties.
#'
#' @importFrom MASS glm.nb polr
#' @importFrom stats constrOptim glm.control rnbinom rnorm
#' @importFrom pscl zeroinfl zeroinfl.control
#'
#' @param X Design matrix (n x q). If NULL, simulated data is generated.
#' @param Y Response vector. If NULL, simulated data is generated.
#' @param beta_ini,phi_ini Initial values; if NULL, estimated from NB or ZINB.
#' @param maxiter Maximum EM iterations (default 200).
#' @param tol Convergence tolerance (default 1e-3).
#' @param ntune Number of tuning parameter combinations (currently internal grid).
#' @return A list with estimated coefficients, selected inflation points,
#'         and variable selection indicators.
#' @export
HMINB <- function(X = NULL, Y = NULL, beta_ini=NULL, phi_ini=NULL, maxiter = 200, tol = 1e-03, ntune = 10){

  if (is.null(X)) {
    set.seed(123)
    n <- 1000
    beta <- c(2, -2, -1.5, 1, 1.5, 0, 0)
    phi <- 1
    K <- c(0, 3, 5, 6, 10, 15, 30)
    gamma1 <- c(3, -5, 4, 6, 0, 0)
    u <- c(0.05, 0.20, 0.30, 0.48, 0.55, 0.64, 0.70)

    o <- produce1(n, beta, phi, K, u, gamma1)
    X <- o$X
    Y <- o$Y
  }

  n <- dim(X)[1]
  X <- cbind(rep(1, n), X)
  q <- dim(X)[2]
  count <- table(Y)
  K.ini <- as.numeric(names(which(count > 2))) #初始计数大于2的视为初始膨胀点///
  M <- length(K.ini)
  #print(table(Y))

  #拟合负二项分布或ZINB分布
  if(is.null(beta_ini)){
    NB <- NULL
    tryCatch({
      NB <- glm.nb(Y~., data = data.frame(X[,-1],Y), control = glm.control(maxit = 500))
      beta.ini <- NB$coefficients
      #print(beta.ini)
      rho2 <- 1.0/abs(beta.ini)
      rho2[1] <- 0
      phi.ini <- max(1.0/NB$theta, 1e-5)
    },error = function(e){print('fitting NB distribution for initialization failed, using ZINB distribution instead')},finally={})
    if (is.null(NB)) {
      suppressWarnings({
        ZINB <- zeroinfl(Y~. | 1, data = data.frame(X[,-1], Y), dist = "negbin", link = "probit", control = zeroinfl.control(method = "CG"))
      })
      beta.ini <- ZINB$coefficients$count
      rho2 <- 1.0/abs(beta.ini)
      rho2[1] <- 0
      phi.ini <- 1/ZINB$theta
    }
  }else{
    beta.ini <- beta_ini
    rho2 <- 1.0/abs(beta.ini)
    rho2[1] <- 0
    phi.ini <- phi_ini
  }


  #拟合累积logit模型
  breaks <- c(-Inf, K.ini, +Inf)
  labels <- 1:(length(K.ini) + 1)
  cut_factor <- cut(Y, breaks = breaks, labels = labels)
  ordered_factor <- factor(cut_factor, ordered = TRUE, levels = labels)

  res <- try(
    polr(
      ordered_factor ~ X[,-1],
      data  = data.frame(X[,-1], ordered_factor),
      Hess  = TRUE
    ),
    silent = TRUE
  )
  if (inherits(res, "try-error")) {
    K.ini <- as.numeric(names(which(count > 2)))
    M <- length(K.ini)
    breaks <- c(-Inf, K.ini, +Inf)
    labels <- 1:(length(K.ini) + 1)
    cut_factor <- cut(Y, breaks = breaks, labels = labels)
    ordered_factor <- factor(cut_factor, ordered = TRUE, levels = labels)
    logit <- polr(ordered_factor~(X[,-1]), data = data.frame(X[,-1], ordered_factor), Hess=TRUE)
    message("Some categories have too few samples; reduce the number of inflation points to merge categories.")
  } else {
    logit <- res
    summary(logit)
  }

  gamma1.ini <- -as.numeric(logit$coefficients)
  rho1 <- 1.0/abs(gamma1.ini)
  gamma0.ini <- as.numeric(logit$zeta)
  u.ini <- log(exp(gamma0.ini)+1)



  #初始化参数
  miu <- 1/(u.ini - c(0,u.ini[-length(u.ini)])) #惩罚系数
  names(gamma1.ini) <- paste0("gamma1$", 1:(q - 1))
  names(beta.ini) <- names(rho2) <- paste0("beta", 0:(q - 1))


  #设置惩罚系数
  varepsilon <- c(1, 3, 5) #u
  #varepsilon <- c(5) #u
  lambda1_set <- c(1, 3, 5) #gamma1
  #lambda1_set <- c(3) #gamma1
  lambda2_set <- c(1, 3, 5) #beta
  #lambda2_set <- c(3) #beta

  l3 <- length(varepsilon)
  l1 <- length(lambda1_set)
  l2 <- length(lambda2_set)
  punish_par <- matrix(0,l3 * l1 * l2,4)
  time <- 1

  # 创建所有参数组合
  param_grid <- expand.grid(varepsilon = varepsilon, lambda1 = lambda1_set, lambda2 = lambda2_set)
  param_grid$time <- 1:(l1 * l2 * l3)
  param_grid <- param_grid[, c("time", "varepsilon", "lambda1", "lambda2")]
  punish_par <- as.matrix(param_grid)

  bic <- matrix(1e20, l3 * l1 * l2, 1)
  result_temp <- list()
  init_set <- list(beta = beta.ini, phi = phi.ini, gamma1 = gamma1.ini, u = u.ini, K = K.ini)
  #print(init_set)
  time <- 1
  for (m3 in 1:l3) {
    vare <- varepsilon[m3]
    for (m1 in 1:l1) {
      lambda1 <- lambda1_set[m1]
      for (m2 in 1:l2) {
        lambda2 <- lambda2_set[m2]
        #print(c(m1,m2,m3))
        #进行第一次优化
        temp <- EM_pen(X, Y, init_set, miu, vare, lambda1, rho1, lambda2, rho2)
        beta <- temp$beta
        phi <- temp$phi
        gamma1 <- temp$gamma1
        u <- temp$u
        K <- temp$K

        result_temp[[time]] <- temp
        #print(temp)
        du <- u - c(0,u[-length(u)])

        #计算似然函数
        lik <- L(X,Y,K,beta,phi,gamma1,u)
        bic[time] <- -2 * lik + (sum(beta!=0) + sum(gamma1!=0) + sum(du!=0)) * log(n)

        time <- time + 1
      }
    }
  }
  #print(result_temp)
  #print(bic)
  id <- which.min(bic)
  output <- result_temp[[id]]
  beta <- output$beta
  phi <- output$phi
  gamma1 <- output$gamma1
  u <- output$u
  K <- output$K

  #选择非零参数
  indicator1 <- which(beta!=0)
  #print(indicator1)
  X1 <- X[,indicator1]
  beta <- beta[indicator1]
  indicator2 <- which(gamma1!=0)
  #print(indicator2)
  X2 <- X[,c(1,(1+indicator2))]
  gamma1 <- gamma1[indicator2]
  set <- list(beta = beta, phi = phi, gamma1 = gamma1, u = u, K = K)
  #cat("dim(X):", dim(X), "\n")
  #cat("dim(X1):", dim(X1), "\n")
  #cat("length(beta):", length(beta), "\n")
  #去掉惩罚项refit
  #print(set)
  cons <- EM_refit(X1, X2, Y, set)
  beta <- cons$beta
  phi <- cons$phi
  gamma1 <- cons$gamma1
  u <- cons$u
  K <- cons$K
  #print(table(Y))
  beta_all <- rep(0,dim(X)[2])
  gamma1_all <- rep(0,dim(X)[2]-1)
  beta_all[indicator1] <- beta
  gamma1_all[indicator2] <- gamma1

  output <- list(indicator1 = (indicator1 - 1)[-1], beta = beta, beta_all=beta_all, phi = phi, kappa = K, u = u, indicator2 = indicator2, gamma1 = gamma1, gamma1_all=gamma1_all)

  return(output)
}


#带惩罚EM优化
EM_pen <- function(X, Y, init_set, miu, vare, lambda1, rho1, lambda2, rho2, maxiter = 50, tol = 1e-03){

  beta <- init_set$beta
  phi <- init_set$phi
  gamma1 <- init_set$gamma1
  u <- init_set$u
  K <- init_set$K

  n <- dim(X)[1]
  q <- dim(X)[2]
  de <- i <- 1


  while (de > tol && i < maxiter && length(u)!=1) {

    temphi <- phi
    tembeta <- beta
    temgamma1 <- gamma1
    temu <- u
    M <- length(u)

    #更新z
    z <- z_upd(X, X, Y, beta, phi, gamma1, u, K)

    #更新beta/phi
    u1 <- matrix(c(1, -1), nrow = 2, ncol = 1)
    c1 <- c(1e-5, -100)

    beta <- beta_upd(X, Y, beta, phi, lambda2, rho2, z)
    #print(beta)
    phi0 <- as.numeric(phi)
    res <- constrOptim(theta = phi0, f = f1, grad = g1, ui = u1, ci = c1, X = X, Y = Y, beta = beta, u = u, z = z, control = list())
    phi <- res$par
    #print(phi)
    #更新gamma
    t <- 1 / norm(X, "2")^2  # 步长
    M <- length(u)
    u <- u_upd(X = X, Y = Y, z = z, gamma1 = gamma1, M = M, vare = vare, miu = miu, u)$par
    #print(u)
    #print(K)
    gamma1 <- gamma1_upd(X, Y, z, K, beta, phi, u, gamma1, t, lambda1, rho1, vare, miu)
    #print(gamma1)

    de <- max(max(abs(beta - tembeta)), max(abs(phi - temphi)), max(abs(u - temu)), max(abs(gamma1 - temgamma1)))
    i <- i + 1
    true <- which((u - c(0,u[-length(u)])) > 0.005)
    K <- K[true]
    u <- u[true]
    miu <- miu[true]
    #print(K)

  }

  beta[abs(beta) < 1e-3] <- 0
  gamma1[abs(gamma1) < 1e-3] <- 0
  output = list(beta = beta, phi = phi, u = u, gamma1 = gamma1, K = K)

  return(output)

}

#不带惩罚优化
EM_refit <- function(X1, X2, Y, init_set, maxiter = 50, tol = 1e-03){

  beta <- init_set$beta
  phi <- init_set$phi
  gamma1 <- init_set$gamma1
  u <- init_set$u
  K <- init_set$K

  vare <- 0
  miu <- rep(0,length(u))
  lambda1 <- 0
  rho1 <- rep(0,length(gamma1))
  lambda2 <- 0
  rho2 <- rep(0,length(beta))
  n <- length(Y)
  M <- length(u)
  de <- i <- 1

  while (de > tol && i < maxiter) {

    temphi <- phi
    tembeta <- beta
    temgamma1 <- gamma1
    temu <- u

    #更新z
    z <- z_upd(X1, X2, Y, beta, phi, gamma1, u, K)

    #更新beta/phi
    u1 <- matrix(c(1, -1), nrow = 2, ncol = 1)
    c1 <- c(1e-5, -100)

    beta <- beta_upd(X1, Y, beta, phi, lambda2, rho2, z)
    #print(beta)
    phi0 <- as.numeric(phi)
    res <- constrOptim(theta = phi0, f = f1, grad = g1, ui = u1, ci = c1, X = X1, Y = Y, beta = beta, u = u, z = z, control = list())
    phi <- res$par

    t <- 1 / norm(X2, "2")^2  # 步长
    #更新gamma
    M <- length(u)
    u <- u_upd(X = X2, Y = Y, z = z, gamma1 = gamma1, M = M, vare = vare, miu = miu, u)$par
    #print(u)
    #print(K)
    gamma1 <- gamma1_upd(X2, Y, z, K, beta, phi, u, gamma1, t, lambda1, rho1, vare, miu)
    #print(gamma1)

    de = max(max(abs(beta - tembeta)), max(abs(phi - temphi)), max(abs(u - temu)), max(abs(gamma1 - temgamma1)))
    i <- i + 1


  }
  #true <- c(1,which((u[-1] - u[1:(length(u) - 1)]) > 0.005) + 1)
  #K <- K[true]
  #u <- u[true]

  output <- list(beta = beta, phi = phi, u = u, gamma1 = gamma1, K=K)
  return(output)

}

#L为似然函数
L <- function(X, Y, K, beta, phi, gamma1, u){

  z <- z_upd(X, X, Y, beta, phi, gamma1, u, K)
  n <- length(Y)
  M <- length(u)

  miu <- exp(X %*% beta)
  p <- matrix(0, n, (M + 1))
  if(length(gamma1) == 0){
    slope <- rep(0,n)
  }else{
    slope <- colSums(t(X[,-1]) * gamma1)
  }
  tem <- outer(exp(-slope), u, FUN = function(a, b) b / (b + a))
  p[, 1] <- tem[,1]
  p[, (M + 1)] <- 1 - tem[,M]
  if(M > 1){
    p[, 2:M] <- tem[, 2:M] - tem[, 1:(M - 1)]
  }
  p[p<=0] <- 0.000001

  a <- lgamma(Y + 1.0/phi) - lgamma(Y + 1) - lgamma(1.0/phi) + Y * log(phi * miu) - (Y + 1.0/phi) * log(1 + phi * miu)

  L <- sum(log(p) * z) + sum(z[,(M + 1)] * a)
  return(L)

}

#z的更新（z指示每一个样本所属类别及概率）
z_upd <- function(X1, X2, Y, beta, phi, gamma1, u, K){
  n <- length(Y)
  M <- length(K)
  #cat("dim(X):", dim(X1), "\n")
  #cat("dim(X1):", dim(X2), "\n")
  #cat("length(beta):", length(beta), "\n")
  #更新delta,p
  X1 <- as.matrix(X1)
  lambda <- as.vector(exp(X1 %*% beta))
  delta <- matrix(0, n, (M + 1))
  for (i in 1:n) {
    delta[i, ifelse(sum(Y[i]==K), which(Y[i]==K), (M + 1))] = 1
  }
  #print(111)
  p <- matrix(0, n, (M + 1))
  if(length(gamma1) == 0){
    slope <- rep(0,n)
  }else{
    slope <- colSums(t(X2[,-1]) * gamma1)
  }
  tem <- outer(exp(-slope), u, FUN = function(a, b) {
    exp_b <- exp(b)
    (exp_b - 1) / (exp_b + a - 1)
  })
  p[, 1] <- tem[,1]
  p[, (M + 1)] <- 1 - tem[,M]
  if(M > 1){
    p[, 2:M] <- tem[, 2:M] - tem[, 1:(M - 1)]
  }
  p[p<=0] <- 0.000001

  #更新z
  z <- matrix(0, n, (M + 1))
  for (m in 1:M) {
    z[,m] <- delta[,m] * p[,m]/(p[,m] + p[,(M+1)] * (exp(lgamma(Y + 1/phi) - lgamma(Y + 1) - lgamma(1/phi) + Y * log(lambda * phi/(1 + lambda * phi)) - log(1 + lambda * phi)/phi)))
  }
  if(M == 1){
    z[,(M + 1)] = 1 - z[,1]
  }else{
    z[,(M + 1)] = 1 - rowSums(z[,1:M])
  }


  return(z)

}


#beta的更新
beta_upd <- function(X, Y, beta, phi, lambda2, rho2, z, maxiter = 200){

  n <- length(Y)
  p <- dim(X)[2]
  z_vec <- z[,ncol(z)]
  beta_old <- beta

  iter_beta <- 0
  repeat{
    abs_sum <- 0
    #print(beta_old)
    for (j in 1:p) {
      mu_vec <- exp(drop(X %*% beta_old) )

      X[,j][X[,j]==0] <- 1e-04#0.5*(1 + 2/pi*atan(X[,j][X[,j]==0]^2/10^-4))
      A_vec <- (Y - mu_vec)/(mu_vec * (1 + mu_vec * phi )) * mu_vec * X[,j]
      B_vec <- (Y - mu_vec)/(mu_vec * (1 + mu_vec * phi )) * mu_vec * X[,j]^2 + (-2 * phi * Y * mu_vec + mu_vec^2 * phi - Y)/(mu_vec^2 * (1 + mu_vec*phi)^2) *mu_vec^2 * X[,j]^2
      #print(sum(z_vec*B_vec))
      omega_vec <- B_vec
      tau_vec <- beta_old[j] - A_vec/B_vec
      #print(sum(A_vec))
      #print(sum(B_vec))

      beta_tilde <- sum( z_vec * (B_vec * beta_old[j] - A_vec))/sum(z_vec * omega_vec)
      #print(beta_tilde)
      dis <- (sum( z_vec * ( - A_vec)) + lambda2 * rho2[j] * sign(beta_tilde))/sum(z_vec * omega_vec)
      betanew <- beta_old[j] + sign(dis) * min(abs(dis),5)
      betanew <- ifelse( (sign(beta_tilde) * (betanew)) >0,  betanew, 0 )

      abs_sum <- abs_sum + abs(beta_old[j] - betanew)
      beta_old[j] <- betanew

    }

    iter_beta <- iter_beta + 1
    if (iter_beta > maxiter) {break}
    if (abs_sum < 10^(-3)) {break}
  }
  return(beta_old)
}


# 软阈值函数
soft_threshold <- function(z, threshold) {
  sign(z) * pmax(abs(z) - threshold, 0)
}
#gamma的更新(近端法)
gamma1_upd <- function(X, Y, z, K, beta, phi, u, gamma1, t, lambda1, rho1, vare, miu, maxiter = 200, tol = 1e-03){

  if(length(gamma1) == 0) return(gamma1)

  n <- length(Y)
  q <- dim(X)[2] - 1
  gamma1_old <- gamma1
  M <- length(u)
  y_k <- gamma1


  for(k in 1:maxiter){
    gamma1_old <- gamma1

    gradient <- g2(X, Y, z, u, gamma1, M, vare, miu)

    y_tem <- y_k - t * gradient

    # 近端操作（带权重的软阈值）
    lambda_t <- lambda1 * t
    gamma1 <- soft_threshold(y_tem, lambda_t * rho1)

    # 更新动量项
    y_k <- gamma1 + (k / (k + 3)) * (gamma1 - gamma1_old)

    # 检查收敛
    if (sum((gamma1 - gamma1_old)^2) < tol^2) break
  }

  return(gamma1_old)

}


#对数障碍函数及其梯度
log_barrier <- function(u, t){

  if (any(u - c(0,u[-length(u)]) < 0)) {
    return(Inf)
  }
  return(-1/t * (sum(log(u - c(0,u[-length(u)])))))

}
grad_log_barrier <- function(u, t){

  grad <- rep(0, length(u))
  if(length(grad)>2){
    grad[1] <- 1/(t * (u[2] - u[1])) - 1/(t * u[1])
    grad[2:(length(u) - 1)] <- 1/t * (1/(u[3:length(u)] - u[2:(length(u) - 1)]) - 1/(u[2:(length(u) - 1)] - u[1:(length(u) - 2)]))
    grad[length(u)] <- -1/(t * (u[length(u)] - u[length(u) - 1]))
  }else if(length(grad)==2){
    grad[1] <- 1/(t * (u[2] - u[1])) - 1/(t * u[1])
    grad[length(u)] <- -1/(t * (u[length(u)] - u[length(u) - 1]))
  }else if(length(grad)==1){
    grad[1] <- - 1/(t * u[1])
  }


  if (any(u - c(0,u[-length(u)]) < 0)) {
    return(rep(Inf, length(u)))
  }
  return(grad)

}
#u的更新（带对数障碍的BFGS拟牛顿方法）
u_upd <- function(X, Y, z, gamma1, M, vare, miu, u, t = 1, tol = 1e-6, max_iter = 100) {

  # 目标函数的组合（目标函数 + 对数障碍项）
  combined_objective_function <- function(x) {
    f_value <- f3(X = X, Y = Y, z = z, u = x, gamma1 = gamma1, M = M, vare = vare, miu = miu)
    barrier_value <- tryCatch(log_barrier(x, t), error = function(e) Inf)
    return(f_value + barrier_value)
  }

  # 目标函数的梯度的组合（目标函数的梯度 + 对数障碍项的梯度）
  combined_gradient_function <- function(x) {
    grad <- g3(X = X, Y = Y, z = z, u = x, gamma1 = gamma1, M = M, vare = vare, miu = miu)
    barrier_grad <- tryCatch(grad_log_barrier(x, t), error = function(e) rep(Inf, length(x)))
    return(grad + barrier_grad)
  }

  initial_point <- u
  # 初始设置
  x <- initial_point
  n <- length(x)
  I <- diag(n)  # 单位矩阵
  H_inv <- I    # Hessian 矩阵的逆的初始值
  alpha <- 1    # 初始步长

  for (iter in 1:max_iter) {

    # 计算组合后的目标函数值和梯度
    combined_value <- combined_objective_function(x)
    combined_grad <- combined_gradient_function(x)

    # 计算搜索方向
    p <- as.numeric(- H_inv %*% as.matrix(combined_grad))

    # 线搜索(二分)
    step_size <- alpha

    while (step_size > 1e-12) {
      new_x <- x + step_size * p
      if (any(!is.finite(new_x))) { step_size <- step_size / 2; next }

      new_value <- combined_objective_function(new_x)
      grad_dot_dir <- sum(combined_grad * p)

      if (is.finite(new_value) && is.finite(grad_dot_dir) &&
          new_value <= combined_value + 0.1 * step_size * grad_dot_dir) {
        break
      }
      step_size <- step_size / 2
    }
    if (step_size <= 1e-12) warning("line search failed")

    # 更新参数
    x_new <- x + step_size * p
    combined_grad_new <- combined_gradient_function(x_new)

    # 近似更新 Hessian 矩阵的逆
    s <- as.numeric(x_new - x)
    y <- as.numeric(combined_grad_new - combined_grad)
    s_mat <- as.matrix(s)
    y_mat <- as.matrix(y)

    ys <- as.numeric(t(y) %*% s)
    eps_pos <- 1e-8
    if (ys > eps_pos) {
      rho <- 1 / ys
      s_mat <- as.matrix(s)
      y_mat <- as.matrix(y)
      H_inv <- (I - rho * (s_mat %*% t(y_mat))) %*% H_inv %*% (I - rho * (y_mat %*% t(s_mat))) + rho * (s_mat %*% t(s_mat))
    } else {
      H_inv <- I
    }

    # 终止条件
    if (sqrt(sum((x_new - x)^2)) < tol) { x <- x_new; break }

    # 更新参数
    x <- x_new
    t <- 2*t

  }

  return(list(par = as.numeric(x), value = combined_objective_function(x)))
}

#phi的优化函数和导函数
f1 <- function(phi, X, Y, beta, u, z){

  y <- Y
  J <- length(u)
  miu <- exp(drop(X %*% beta))
  a <- lgamma(y + 1.0/phi) - lgamma(y + 1) - lgamma(1.0/phi) + y * log(phi * miu) - (y + 1.0/phi) * log(1 + phi * miu)
  obj <- -1.0*sum(z[,J+1]*a)

  return(obj)
}

g1 <- function(phi, X, Y, beta, u, z){

  y <- Y
  J <- length(u)
  miu <- exp(drop(X%*%beta))
  g <- digamma(y + 1.0/phi) * (-1.0/phi^2) - digamma(1.0/phi) * (-1.0/phi^2) + y/phi + 1.0/phi^2 * log(1 + miu * phi)-(y + 1.0/phi) * (miu/(1 + miu*phi))
  g <- -sum(z[,J+1] * g)

  return(c(g))
}

#gamma的优化函数和导数
f2 <- function(X, Y, z, u, gamma1, M, vare, miu){

  n=dim(X)[1]

  p=matrix(0,n,(M+1))
  if(length(gamma1) == 0){
    slope <- rep(0,n)
  }else{
    slope <- colSums(t(X[,-1]) * gamma1)
  }
  tem <- outer(exp(-slope), u, FUN = function(a, b) {
    exp_b <- exp(b)
    (exp_b - 1) / (exp_b + a - 1)
  })
  p[, 1] <- tem[,1]
  p[, (M + 1)] <- 1 - tem[,M]
  if(M > 1){
    p[, 2:M] <- tem[, 2:M] - tem[, 1:(M - 1)]
  }
  p[p<=0]=0.000001

  f=-sum(log(p)*z)+vare*sum(miu*(u-c(0,u[-length(u)])))

  return(f)

}

g2 <- function(X, Y, z, u, gamma1, M, vare, miu){

  n=dim(X)[1]
  q=length(gamma1)

  p=matrix(0,n,(M+1))
  if(length(gamma1) == 0){
    slope <- rep(0,n)
  }else{
    slope <- colSums(t(X[,-1]) * gamma1)
  }
  tem <- outer(exp(-slope), u, FUN = function(a, b) {
    exp_b <- exp(b)
    (exp_b - 1) / (exp_b + a - 1)
  })
  p[, 1] <- tem[,1]
  p[, (M + 1)] <- 1 - tem[,M]
  if(M > 1){
    p[, 2:M] <- tem[, 2:M] - tem[, 1:(M - 1)]
  }
  p[p<=0]=0.000001

  # ---------- 对 gamma1 的梯度 ----------
  if (q == 0) return(numeric(0))

  tem = outer(exp(-slope), u, FUN = function(a, b) {
    c = exp(b) - 1
    d = c + a
    a * c / (d^2)
  })

  tem2 = matrix(0, n, M + 1)
  tem2[, 1] = tem[, 1]
  tem2[, M + 1] = -tem[, M]
  if (M > 1) {
    tem2[, 2:M] = tem[, 2:M] - tem[, 1:(M-1)]
  }
  T_vec = rowSums((z / p) * tem2)

  g2 = rep(0,q)
  for (k in 2:(q+1)) {
    g2[k-1] = -sum(X[, k] * T_vec)
  }

  return(g2)
}

#u的优化函数和导函数
f3 <- function(X, Y, z, u, gamma1, M, vare, miu){
  u <- pmin(u, 30)
  n <- length(Y)
  p <- matrix(0, n, (M + 1))
  if(length(gamma1) == 0){
    slope <- rep(0,n)
  }else{
    slope <- colSums(t(X[,-1]) * gamma1)
  }
  tem <- outer(exp(-slope), u, FUN = function(a, b) {
    exp_b <- exp(b)
    (exp_b - 1) / (exp_b + a - 1)
  })
  p[, 1] <- tem[,1]
  p[, (M + 1)] <- 1 - tem[,M]
  if(M > 1){
    p[, 2:M] <- tem[, 2:M] - tem[, 1:(M - 1)]
  }
  p[p<=0] <- 0.000001

  f <- -sum(log(p) * z) + vare * sum(miu * (u - c(0,u[-length(u)])))

  return(f)

}

g3 <- function(X, Y, z, u, gamma1, M, vare, miu){
  u <- pmin(u, 30)
  n <- length(Y)
  q <- length(gamma1)
  p <- matrix(0, n, (M + 1))
  if(length(gamma1) == 0){
    slope <- rep(0,n)
  }else{
    slope <- colSums(t(X[,-1]) * gamma1)
  }
  tem <- outer(exp(-slope), u, FUN = function(a, b) {
    exp_b <- exp(b)
    (exp_b - 1) / (exp_b + a - 1)
  })
  p[, 1] <- tem[,1]
  p[, (M + 1)] <- 1 - tem[,M]
  if(M > 1){
    p[, 2:M] <- tem[, 2:M] - tem[, 1:(M - 1)]
  }
  p[p<=0] <- 0.000001

  tem2 <- outer(exp(-slope), u, FUN = function(a, b) {
    exp_b <- exp(b)
    a * exp_b / (exp_b + a - 1)^2
  })
  tem2 <- tem2 * (z[,-1]/p[,-1] - z[,(1:M)]/p[,(1:M)])
  g <- colSums(tem2)
  pen <- rep(0, M)
  if (M == 1) {
    pen[1] <- vare * miu[1]
  } else {
    pen[1:(M-1)] <- vare * (miu[1:(M-1)] - miu[2:M])
    pen[M] <- vare * miu[M]
  }
  g <- g+pen

  return(g)

}


produce1 <- function(n, beta, phi, K, u, gamma1){

  q <- length(beta) - 1
  M <- length(K)
  X1 <- matrix(data = rnorm((n/2) * q, mean = 0, sd = 2), nrow = n/2, ncol = q)
  X2 <- matrix(data = rnorm((n/2) * q, mean = 2, sd = 2), nrow = n/2, ncol = q)
  X <- rbind(X1, X2)
  Y <- rep(0, n)
  p <- rep(0, M + 1) # 注意：膨胀概率有M+1个

  siz <- 1/phi
  for (i in 1:n) {
    pro <- 1/(1 + exp(beta[1] + sum(X[i,] * beta[-1])) * phi)
    slope <- sum(gamma1 * X[i,])
    tem <- (exp(u) - 1) / (exp(u) + exp(-slope) - 1)

    p[1] <- tem[1]
    p[M + 1] <- 1 - tem[M]
    if(M > 1){
      p[2:M] <- tem[2:M] - tem[1:(M - 1)]
    }

    # 确保概率非负且和为1
    p <- pmax(p, 0)
    p <- p / sum(p)

    # 生成负二项分布的计数
    y1 <- rnbinom(1, size = siz, prob = pro)

    # 从混合分布中抽样
    Y[i] <- sample(c(K, y1), size = 1, replace = TRUE, prob = p)
  }

  return(list(X = X, Y = Y))
}
