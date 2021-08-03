Gibbs_sampler <- function(
  T_obs, # Observed recurrent event times
  n, # Vector with the number of recurrent events of each individual
  c_i, # Minimum of survival and censoring time
  cens = rep(FALSE, nrow(T_obs)), # Censoring indicator
  x = matrix(NA_real_, nrow = nrow(T_obs), ncol = 0), # Covariates
  n_iter = 1e2L, # Number of iterations of the Gibbs sampler
  burnin = as.integer(n_iter / 4), # Number of burnin iterations
  thin = 1L, # Thinning rate
  a_lambda = 1, # Prior parameters for `lambda`
  b_lambda = 1
) {
  
  if (!is.matrix(x)) stop("This code requires x to be a matrix.")
  
  # Number of recorded iterations
  n_recorded <- as.integer((n_iter - burnin) / thin)
  
  # Number of individuals
  L <- length(n)
  
  Y_obs <- log(t(diff(t(cbind(0, T_obs)))))
  n_max <- max(n)
  
  # Number of covariates
  q <- dim(x)[2]
  if (q == 1) stop("This code does not work with exactly one covariate.")
  
  # The parameter "m" in Neal's Algorithm 8
  m_Neal8 <- 2L
  
  
  ## Specify hyperparameters.
  sigma2_beta <- 10^2
  sigma2_gamma <- 10^2
  
  sigma2_m <- 100
  sigma2_delta <- 100
  
  # Gamma(a_M, b_M) prior on the DP concentration parameter M.
  a_M <- 2
  b_M <- 1
  
  a_r <- 1
  b_r <- 1
  
  nu_sigma2 <- 4.02
  sigma2_0 <- 2.02 / 4.02
  
  nu_eta2 <- 4.02
  eta2_0 <- 2.02 / 4.02
  
  a_eta2 <- (nu_eta2 + L) / 2
  
  
  ## Objects to store the Gibbs samples
  
  # Regression coefficients
  beta_Gibbs <- matrix(NA_real_, n_recorded, q)
  gamma_Gibbs <- matrix(NA_real_, n_recorded, q)
  
  # Number of recurrent events
  N_Gibbs <- matrix(NA_integer_, n_recorded, L)
  
  # Unobserved log gap times
  # We will increase the dimension of this array if the size of N demands.
  Y_unobs_Gibbs <- array(NA_real_, dim = c(n_recorded, sum(cens), 0))
  
  # Survival times
  S_Gibbs <- matrix(NA_real_, n_recorded, L)
  
  # Cluster allocation
  s_Gibbs <- matrix(NA_integer_, n_recorded, L)
  
  # Number of clusters
  K_Gibbs <- rep(NA_integer_, n_recorded)
  
  # Cluster specific parameters
  m_Gibbs <- array(NA_real_, dim = c(n_recorded, L, 2))
  delta_Gibbs <- matrix(NA_real_, n_recorded, L)
  
  # DP concentration parameter M
  M_Gibbs <- rep(NA_real_, n_recorded)
  
  # Variance parameters
  sigma2_Gibbs <- rep(NA_real_, n_recorded)
  eta2_Gibbs <- rep(NA_real_, n_recorded)
  r_Gibbs <- rep(NA_real_, n_recorded)
  lambda_Gibbs <- rep(NA_real_, n_recorded)
  
  
  ## Initialize the Gibbs sampler.
  N <- n
  N_max <- n_max
  
  # Set lambda to the average of the observed N,
  # or, if there aren't any, then the average of n.
  lambda <- mean(n[if (all(cens)) 1:L else !cens])
  
  # Compute the initial value for r to match the variance of the observations.
  r <- max(
    0.01,
    lambda^2 / (var(n[if (sum(!cens) < 2) 1:L else !cens]) - lambda)
  )
  
  S <- ifelse(
    test = cens,
    yes = 1.5,
    no = 1
  ) * c_i
  
  M <- 1
  
  K <- 1L # Number of clusters
  s <- rep(1L, L) # Cluster allocation
  s_n <- as.vector(table(s))
  
  # Cluster specific parameters
  m <- matrix(0, L, 2)
  delta <- rep(NA_real_, L)
  
  for (h in 1:K) {
    index <- which(s == h)
    m[index, 1] <- mean(Y_obs[index, ], na.rm = TRUE)
    delta[index] <- mean(log(S[index]))
  }
  
  # Regression coefficients
  beta <- rep(0, q)
  gamma <- rep(0, q)
  
  # Variances
  sigma2 <- var(as.vector(Y_obs), na.rm = TRUE) / 2
  eta2 <- var(log(S)) / 2
  if (eta2 == 0) eta2 <- eta2_0
  
  Y_unobs <- matrix(NA_real_, nrow = L, ncol = 0)
  Y <- Y_obs
  
  T_n <- numeric(L)
  for (i in 1:L) if (n[i] > 0L) T_n[i] <- T_obs[i, n[i]]
  T_N <- T_n
  
  
  # Function that computes Sigma_Y divided by sigma2
  Sigma_Y_div_sigma2_compute <- function(m_2, N) {
    
    L <- length(N)
    if (length(m_2) == 1L) m_2 <- rep(m_2, L)
    
    Sigma_Y_div_sigma2 <- list()
    
    for (i in 1:L) {
      
      j_1 <- .row(dim = c(N[i], N[i]))
      j_2 <- .col(dim = c(N[i], N[i]))
      
      Sigma_Y_div_sigma2[[i]] <- (
        m_2[i]^abs(j_1 - j_2) - m_2[i]^(j_1 + j_2)
      ) / (1 - m_2[i]^2)
      
    }
    
    Sigma_Y_div_sigma2
    
  }
  
  Sigma_Y_div_sigma2 <- Sigma_Y_div_sigma2_compute(m[, 2], N)
  
  
  # Function that computes log(sum(exp(vec)))
  # while avoiding over- and underflow errors
  # We rename this function to avoid spending time on repeatedly exporting
  # this function from the package `matrixStats`.
  log_sum_exp <- matrixStats::logSumExp
  
  # Function that computes LS_1 and LS_2
  LS_compute <- function(Sigma_Y_div_sigma2, sigma2) {
    
    L <- length(Sigma_Y_div_sigma2)
    LS <- matrix(NA_real_, nrow = L, ncol = 2)
    
    for(i in 1:L) {
      
      N_i <- nrow(Sigma_Y_div_sigma2[[i]])
      
      if (N_i > 0) {
        
        Sigma_Y_i <- Sigma_Y_div_sigma2[[i]] * sigma2
        tmp <- diag(Sigma_Y_i) / 2
        
        LS[i, 1] <- log_sum_exp(tmp)
        LS[i, 2] <- log_sum_exp(tmp + rep(tmp, rep.int(N_i, N_i)) + Sigma_Y_i)
        
      }
    }
    
    LS
    
  }
  
  
  # The Fenton-Wilkinson approximation to the normalization constants
  # Result is on a log scale.
  # The result is summed over the individuals considered.
  Pr_sum_log <- function(beta, gamma, m_1, delta, eta2, LS, x) {
    
    if (is.vector(LS)) LS <- t(as.matrix(LS))
    
    sum(ifelse(test = is.na(LS[, 1]), yes = 0, no = pnorm(
      q = 0,
      mean = x %*% beta + m_1 + 2 * LS[, 1] - LS[, 2] / 2 - x %*% gamma - delta,
      sd = sqrt(LS[, 2] - 2 * LS[, 1] + eta2),
      log.p = TRUE
    )))
    
  }
  
  
  # Function that samples from a left-truncated normal distribution
  # using the inverse transformation method.
  rnorm_truncated <- function(
    n, mean = 0,
    sd = 1,
    lower_bound = -Inf,
    upper_bound = Inf
  ) {
    
    if (lower_bound == -Inf) return(qnorm(p = pnorm(
      q = upper_bound,
      mean = mean,
      sd = sd,
      lower.tail = TRUE,
      log.p = TRUE
    ) - rexp(n), mean = mean, sd = sd, lower.tail = TRUE, log.p = TRUE))
    
    if (upper_bound == Inf) return(qnorm(p = pnorm(
      q = lower_bound,
      mean = mean,
      sd = sd,
      lower.tail = FALSE,
      log.p = TRUE
    ) - rexp(n), mean = mean, sd = sd, lower.tail = FALSE, log.p = TRUE))
    
    qnorm(p = runif(
      n,
      min = pnorm(q = lower_bound, mean = mean, sd = sd),
      max = pnorm(q = upper_bound, mean = mean, sd = sd)
    ), mean = mean, sd = sd)
    
  }
  
  
  # Function for univariate slice sampling
  # This is an edited version of Radford Neal's code avaialable at
  # (https://www.cs.toronto.edu/~radford/ftp/slice-R-prog)
  slice_sampling <- function(x0, g, w = 1, lower = -Inf, upper = +Inf) {
    
    # UNIVARIATE SLICE SAMPLING WITH STEPPING OUT AND SHRINKAGE
    #
    # Performs a slice sampling update from an initial point to a new point that
    # leaves invariant the distribution with the specified log density function.
    #
    # Arguments:
    #   x0    Initial point
    #   g     Function returning the log probability density plus constant
    #   w     Size of the steps for creating interval
    #   lower Lower bound on support of the distribution
    #   upper Upper bound on support of the distribution
    #
    # The log density function may return -Inf for points outside the support 
    # of the distribution.  If a lower and/or upper bound is specified for the
    # support, the log density function will not be called outside such limits.
    #
    # The value of this function is the new point sampled.
    
    return_x0 <- function() {
      
      warning(
        "Likelihood return NaN while slice sampling. Returning starting point."
      )
      
      return(x0)
      
    }
    
    # Determine the slice level, in log terms.
    logy <- g(x0) - rexp(1)
    if (is.na(logy)) return(return_x0())
    
    # Find the initial interval to sample from.
    u <- runif(1, 0, w)
    
    L <- x0 - u
    
    # The parentheses should guarantee that x0 is in [L,R], even with roundoff:
    R <- x0 + (w - u)
    
    # Expand the interval until its ends are outside the slice.
    repeat {
      if (L <= lower) break
      g_L <- g(L)
      if (is.na(g_L)) return(return_x0())
      if (g_L <= logy) break
      L <- L - w
    }
    
    repeat {
      if (R >= upper) break
      g_R <- g(R)
      if (is.na(g_R)) return(return_x0())
      if (g_R <= logy) break
      R <- R + w
    }
    
    # Shrink interval to lower and upper bounds.
    L <- max(L, lower)
    R <- min(R, upper)
    
    # Sample from the interval, shrinking it on each rejection.
    repeat {
      x1 <- runif(1, L, R)
      g_x1 <- g(x1)
      if (is.na(g_x1)) return(return_x0())
      if (g_x1 >= logy) return(x1)
      if (x1 > x0) R <- x1 else L <- x1
    }
    
  }
  
  
  pb <- txtProgressBar(max = n_iter, style = 3)
  
  ## Main loop with the Gibbs sampler iterations
  for (iter in 1:n_iter) {
    
    # This follows the steps in Algorithm S1 of the supplemental material.
    
    
    # Step 1: Update beta and gamma.
    LS <- LS_compute(Sigma_Y_div_sigma2, sigma2)
    
    if (q != 0L) {
      
      # We compute `Sigma_beta_star` divided by sigma2.
      Sigma_beta_star_inv <- sigma2 / sigma2_beta * diag(q)
      
      for (i in 1:L) Sigma_beta_star_inv <- Sigma_beta_star_inv + ifelse(
        test = N[i] == 0L, yes = 0, no = 1 + (N[i] - 1) * (1 - m[i, 2])^2
      ) * tcrossprod(x[i, ])
      
      Sigma_beta_star <- solve(Sigma_beta_star_inv)
      mu_beta_star <- numeric(q)
      
      if (N[i] > 0) mu_beta_star <- Sigma_beta_star %*% colSums(
          x * rowSums(cbind(
            Y[, 1] - m[, 1],
            (1 - m[, 2]) * (Y[, -1] - m[, 1] - m[, 2] * (Y[, -N_max] - m[, 1]))
        ), na.rm = TRUE))
      
      # We compute `Sigma_eta_star` divided by eta2.
      Sigma_gamma_star <- solve(eta2 / sigma2_gamma * diag(q) + crossprod(x))
      
      mu_gamma_star <- Sigma_gamma_star %*% crossprod(x, log(S) - delta)
      
    }
    
    for (k in seq_len(q)) {
      
      # Update beta_k.
      tmp <- solve(Sigma_beta_star[-k, -k]) %*% Sigma_beta_star[k, -k]
      
      mu_beta_k_star <- mu_beta_star[k] + sum(
        tmp * (beta[-k] - mu_beta_star[-k])
      )
      
      Sigma_beta_k_star <- sigma2 * (Sigma_beta_star[k, k] - sum(
        tmp * Sigma_beta_star[k, -k]
      ))
      
      beta[k] <- slice_sampling(
        x0 = beta[k],
        g = function(beta_k) {
          
          beta[k] <- beta_k
          
          dnorm(
            x = beta_k,
            mean = mu_beta_k_star,
            sd = sqrt(Sigma_beta_k_star),
            log = TRUE
          ) - Pr_sum_log(beta, gamma, m[, 1], delta, eta2, LS, x)
          
        },
        w = sqrt(Sigma_beta_k_star)
      )
      
      # Update gamma_k.
      tmp <- solve(Sigma_gamma_star[-k, -k]) %*% Sigma_gamma_star[k, -k]
      
      mu_gamma_k_star <- mu_gamma_star[k] + sum(
        tmp * (gamma[-k] - mu_gamma_star[-k])
      )
      
      Sigma_gamma_k_star <- eta2 * (Sigma_gamma_star[k, k] - sum(
        tmp * Sigma_gamma_star[k, -k]
      ))
      
      gamma[k] <- slice_sampling(
        x0 = gamma[k],
        g = function(gamma_k) {
          
          gamma[k] <- gamma_k
          
          dnorm(
            x = gamma_k,
            mean = mu_gamma_k_star,
            sd = sqrt(Sigma_gamma_k_star),
            log = TRUE
          ) - Pr_sum_log(beta, gamma, m[, 1], delta, eta2, LS, x)
          
        },
        w = sqrt(Sigma_gamma_k_star)
      )
      
    }
    
    
    # Step 2: Update N_i, Y_i and S_i for each censored individual.
    for (i in which(cens)) {
      
      
      # Step 2(a): Update N_i using a reversible jump sampler
      
      # Propose N_i from a negative binomial distribution left-truncated at n_i.
      N_proposal <- as.integer(qnbinom(
        p = runif(n = 1, min = pnbinom(q = n[i] - 1, size = r, mu = lambda)),
        size = r,
        mu = lambda
      ))
      
      N_d_proposal <- N_proposal - n[i]
      
      # Propose the unobserved Y.
      T_unobs_i_proposal <- c(T_n[i], rep(NA_real_, N_d_proposal))
      
      for (j in seq_len(N_d_proposal)) T_unobs_i_proposal[j + 1] <- runif(
        n = 1,
        min = if (j == 1) c_i[i] else T_unobs_i_proposal[j],
        max = S[i]
      )
      
      Y_unobs_i_proposal <- log(diff(T_unobs_i_proposal))
      
      mu <- sum(x[i, ] * beta) + m[i, 1]
      
      # Compute the acceptance probability.
      C_div_f_log <- function(Y_unobs_i) {
        
        if (any(is.infinite(Y_unobs_i))) return(-Inf)
        N_d <- length(Y_unobs_i)
        
        LS_tmp <- LS_compute(
          Sigma_Y_div_sigma2_compute(m[i, 2], n[i] + N_d),
          sigma2
        )
        
        Pr <- Pr_sum_log(beta, gamma, m[i, 1], delta[i], eta2, LS_tmp, x[i, ])
        if (N_d == 0L) return(-Pr)
        sigma <- sqrt(sigma2)
        
        if (n[i] == 0L) {
          C_log <- dnorm(
            x = Y_unobs_i[1],
            mean = mu,
            sd = sigma,
            log = TRUE
          ) + sum(dnorm(
            x = Y_unobs_i[-1],
            mean = mu + m[i, 2] * (Y_unobs_i[-N_d] - mu),
            sd = sigma,
            log = TRUE
          )) - Pr
        } else {
          C_log <- sum(dnorm(
            x = Y_unobs_i,
            mean = mu + m[i, 2] * (c(Y_obs[i, n[i]], Y_unobs_i[-N_d]) - mu),
            sd = sigma,
            log = TRUE
          )) - Pr
        }
        
        f_log <- log(
          exp(Y_unobs_i[1]) + T_n[i] - c_i[i]
        ) - log(S[i] - c_i[i]) + sum(
          Y_unobs_i[-1] - log(
            # We use pmax as the following expression can be negative
            # due to floating-point error.
            pmax(S[i] - T_n[i] - cumsum(exp(Y_unobs_i[-N_d])), 0)
          )
        )
        
        C_log - f_log
        
      }
      
      N_d <- N[i] - n[i]
      Y_unobs_i <- Y_unobs[i, seq_len(N_d)]
      a <- exp(C_div_f_log(Y_unobs_i_proposal) - C_div_f_log(Y_unobs_i))
      
      if (is.na(a)) {
        warning("`a` evaluated as NA. Skipping this observation.")
        next
      } else if (runif(1) < a) {
        
        N[i] <- N_proposal
        N_d <- N_d_proposal
        Y_unobs_i <- Y_unobs_i_proposal
        
        Sigma_Y_div_sigma2[[i]] <- Sigma_Y_div_sigma2_compute(
          m[i, 2],
          N_proposal
        )
        
        LS[i, ] <- LS_compute(
          Sigma_Y_div_sigma2_compute(m[i, 2], N_proposal),
          sigma2
        )
        
      }
      
      
      # Step 2(b): Update Y_i using the inverse transformation method.
      
      ## We update Y_unobs_i according to the marginal full conditionals
      ## of its elements,
      ## regardless of whether we accept or reject the proposed (N, Y_unobs).
      for (j in seq_len(N_d)) {
        
        first <- j == 1L
        last <- j == N_d
        
        Y_previous_mu <- ifelse(first, Y_obs[i, n[i]], Y_unobs_i[j - 1L]) - mu
        
        Y_unobs_i[j] <- rnorm_truncated(
          n = 1,
          mean = mu + ifelse(
            test = last,
            yes = ifelse(N[i] == 1, 0, m[i, 2] * Y_previous_mu),
            no = m[i, 2] / (1 + m[i, 2]^2) * (Y_unobs_i[j + 1L] - mu + ifelse(
              test = first,
              yes = 0,
              no = Y_previous_mu
            ))
          ),
          sd = sqrt(ifelse(last, sigma2, sigma2 / (1 + m[i, 2]^2))),
          lower_bound = ifelse(first, log(c_i[i] - T_n[i]), -Inf),
          upper_bound = log(S[i] - T_n[i] - sum(exp(Y_unobs_i[-j])))
        )
        
      }
      
      T_N[i] <- T_n[i] + sum(exp(Y_unobs_i))
      
      if (N_d > ncol(Y_unobs)) {
        
        Y_unobs <- cbind(
          Y_unobs,
          matrix(NA_real_, nrow = L, ncol = N_d - ncol(Y_unobs))
        )
        
        Y_unobs_Gibbs <- array(
          data = c(Y_unobs_Gibbs, rep(NA_real_, n_recorded * sum(cens))),
          dim = c(n_recorded, sum(cens), ncol(Y_unobs))
        )
        
      }
      
      Y_unobs[i, ] <- c(Y_unobs_i, rep(NA_real_, ncol(Y_unobs) - N_d))
      
      
      # Step 2(c): Update the unobserved survival time.
      # The exponentiation can lead to `Inf`.
      # Floating-point error might violate the lower bound.
      # We therefore use max/min.
      tmp <- max(T_N[i], c_i[i])
      
      S[i] <- max(min(exp(rnorm_truncated(
        n = 1,
        mean = sum(x[i, ] * gamma) + delta[i],
        sd = sqrt(eta2),
        lower_bound = log(tmp)
      )), 1e100), tmp)
      
    }
    
    N_max <- max(N)
    Y <- cbind(Y_obs, matrix(NA_real_, nrow = L, ncol = N_max - n_max))
    
    for (i in which(cens)) {
      N_d <- N[i] - n[i]
      Y[i, n[i] + seq_len(N_d)] <- Y_unobs[i, seq_len(N_d)]
    }
    
    Z <- Y - as.vector(x %*% beta)
    
    
    # Step 3: Update m and delta via Neal's Algorithm 8.
    for (i in 1:L) {
      
      ## Update the cluster allocations
      
      # Sample candidate parameter values for new clusters.
      m_proposal <- matrix(rnorm(
        n = 2L * m_Neal8,
        sd = sqrt(sigma2_m)
      ), ncol = 2)
      
      delta_proposal <- rnorm(n = m_Neal8, sd = sqrt(sigma2_delta))
      
      # Remove individual i from their cluster.
      s_i <- s[i]
      s_n[s_i] <- s_n[s_i] - 1L
      
      if (s_n[s_i] == 0L) {
        
        K <- K - 1L # Number of clusters
        
        # Ensure that the clusters are labeled consecutively.
        s[s > s_i] <- s[s > s_i] - 1L
        s_n <- s_n[-s_i]
        
        m_proposal[1, ] <- m[i, ]
        delta_proposal[1] <- delta[i]
        
        # Ensure that `index` defined next identifies the parameters of cluster
        # s_i correctly.
        s[i] <- NA
        
      }
      
      # Keep track of m_i_2 to see whether we need to recompute LS.
      m_i_2 <- m[i, 2]
      
      Z_i <- Z[i, seq_len(N[i])]
      prob_log <- rep(NA_real_, K + m_Neal8)
      index <- match(1:K, s)
      
      # Compute the probability that individual i belongs to cluster h.
      prob_log_compute <- function(m, delta) sum(dnorm(
        x = Z_i,
        mean = c(m[1], m[1] + m[2] * (Z_i[-N[i]] - m[1])),
        sd = sqrt(sigma2),
        log = TRUE
      )) - Pr_sum_log(
        beta = beta,
        gamma = gamma,
        m_1 = m[1],
        delta = delta,
        eta2,
        LS = LS_compute(Sigma_Y_div_sigma2_compute(m[2], N[i]), sigma2),
        x = x[i, ]
      ) + dnorm(
        x = log(S[i]),
        mean = sum(x[i, ] * gamma) + delta,
        sd = sqrt(eta2),
        log = TRUE
      )
      
      for (h in 1:K) {
        prob_log[h] <- prob_log_compute(m[index[h], ], delta[index[h]])
      }
      
      for (h in 1:m_Neal8) {
        prob_log[K + h] <- prob_log_compute(m_proposal[h, ], delta_proposal[h])
      }
      
      prob_log <- log(c(s_n, rep(M / m_Neal8, m_Neal8))) + prob_log
      prob_log_max <- max(prob_log)
      
      s[i] <- sample.int(
        n = K + m_Neal8,
        size = 1,
        prob = if(is.finite(prob_log_max)) {
          exp(prob_log - prob_log_max)
        } else if (prob_log_max == Inf) prob_log == Inf
      )
      
      
      if (s[i] > K) {
        
        m[i, ] <- m_proposal[s[i] - K, ]
        delta[i] <- delta_proposal[s[i] - K]
        
        K <- K + 1L # Number of clusters
        s[i] <- K
        s_n <- c(s_n, 1L)
        
      } else {
        
        s_n[s[i]] <- s_n[s[i]] + 1L
        
        m[i, ] <- m[index[s[i]], ]
        delta[i] <- delta[index[s[i]]]
        
      }
      
      
      # Update LS if m[i, 2] changed.
      if (m[i, 2] != m_i_2) LS[i, ] <- LS_compute(
        Sigma_Y_div_sigma2 = Sigma_Y_div_sigma2_compute(m[i, 2], N[i]),
        sigma2 = sigma2
      )
      
    }
    
    
    # Update all cluster specific parameters.
    tmp <- rowSums(
      cbind(Z[, 1], (1 - m[, 2]) * (Z[, -1] - m[, 2] * Z[, -N_max])),
      na.rm = TRUE
    )
    
    for (h in 1:K) {
      
      index <- which(s == h)
      
      # m_1
      Sigma_m_1 <- 1 / (1 / sigma2_m + sum(
        ifelse(N[index] == 0L, 0, 1 + (N[index] - 1) * (1 - m[index, 2])^2)
      ) / sigma2)
      
      m[index, 1] <- slice_sampling(
        x0 = m[index[1], 1],
        g = function(m_star_1) dnorm(
          x = m_star_1,
          mean = Sigma_m_1 * sum(tmp[index]) / sigma2,
          sd = sqrt(Sigma_m_1),
          log = TRUE
        ) - Pr_sum_log(
          beta, gamma, m_star_1, delta[index], eta2, LS[index, ], x[index, ]
        ),
        w = sqrt(Sigma_m_1)
      )
      
      
      # m_2
      mu_m_2_sum <- 0
      Sigma_m_2_sum <- 0
      
      for (i in index) {
        Z_i <- Z[i, seq_len(N[i])] - m[i, 1]
        mu_m_2_sum <- mu_m_2_sum + sum(Z_i[-1] * Z_i[-N[i]])
        Sigma_m_2_sum <- Sigma_m_2_sum + sum(Z_i[-N[i]]^2)
      }
      
      Sigma_m_2 <- 1 / (1 / sigma2_m + Sigma_m_2_sum / sigma2)
      
      m[index, 2] <- slice_sampling(
        x0 = m[index[1], 2],
        g = function(m_star_2) dnorm(
          x = m_star_2,
          mean = Sigma_m_2 * mu_m_2_sum / sigma2,
          sd = sqrt(Sigma_m_2),
          log = TRUE
        ) - Pr_sum_log(
          beta = beta,
          gamma = gamma,
          m_1 = m[index, 1],
          delta = delta[index],
          eta2,
          LS = LS_compute(
            Sigma_Y_div_sigma2 = Sigma_Y_div_sigma2_compute(m_star_2, N[index]),
            sigma2 = sigma2
          ),
          x = x[index, ]
        ),
        w = sqrt(Sigma_m_2)
      )
      
      
      # delta
      Sigma_delta <- 1 / (1 / sigma2_delta + length(index) / eta2)
      xTgamma <- x[index, ] %*% gamma
      
      delta[index] <- slice_sampling(
        x0 = delta[index[1]],
        g = function(delta_star) dnorm(
          x = delta_star,
          mean = Sigma_delta * sum(log(S[index]) - xTgamma ) / eta2,
          sd = sqrt(Sigma_delta),
          log = TRUE
        ) - Pr_sum_log(
          beta, gamma, m[index, 1], delta_star, eta2, LS[index, ], x[index, ]
        ),
        w = sqrt(Sigma_delta)
      )
      
    }
    
    
    # Step 4: Update the concentration parameter M
    # following Escobar & West (1995, Equation 13).
    eta_M <- rbeta(1, M + 1, L)
    
    M <- rgamma(
      n = 1,
      shape = a_M + K - (
        runif(1) < 1 / (1 + (a_M + K - 1) / (L * (b_M - log(eta_M))))
      ),
      rate = b_M - log(eta_M)
    )
    
    
    # Step 5: Update sigma2, eta2 and lambda.
    
    # Update sigma2.
    tmp <- matrix(NA_real_, nrow = L, ncol = N_max)
    tmp[, 1] <- Z[, 1] - m[, 1]
    tmp[, -1] <- Z[, -1] - m[, 1] - m[, 2] * (Z[, -N_max] - m[, 1])
    
    a_sigma2 <- (nu_sigma2 + sum(N)) / 2
    b_sigma2 <- (nu_sigma2 * sigma2_0 + sum(tmp * tmp, na.rm = TRUE)) / 2
    
    Sigma_Y_div_sigma2 <- Sigma_Y_div_sigma2_compute(m[, 2], N)
    
    sigma2 <- 1 / slice_sampling(
      x0 = 1 / sigma2,
      g = function(sigma2_inv) dgamma(
        x = sigma2_inv,
        shape = a_sigma2,
        rate = b_sigma2,
        log = TRUE
      ) - Pr_sum_log(
        beta = beta,
        gamma = gamma,
        m_1 = m[, 1],
        delta = delta,
        eta2,
        LS = LS_compute(Sigma_Y_div_sigma2, 1 / sigma2_inv),
        x = x
      ),
      w = sqrt(a_sigma2) / b_sigma2,
      lower = 0
    )
    
    
    # Update eta2.
    mu <- x %*% gamma + delta
    
    b_eta2 <- (nu_eta2 * eta2_0 + sum((log(S) - mu)^2)) / 2
    
    eta2 <- 1 / slice_sampling(
      x0 = 1 / eta2,
      g = function(eta2_inv) dgamma(
        x = eta2_inv,
        shape = a_eta2,
        rate = b_eta2,
        log = TRUE
      ) - Pr_sum_log(beta, gamma, m[, 1], delta, 1 / eta2_inv, LS, x),
      w = sqrt(a_eta2) / b_eta2,
      lower = 0
    )
    
    
    # Update r and lambda.
    r <- slice_sampling(
      x0 = r,
      g = function(r) dgamma(
        x = r, shape = a_r, rate = b_r, log = TRUE
      ) + sum(dnbinom(x = N, size = r, mu = lambda, log = TRUE)),
      w = sqrt(a_r) / b_r,
      lower = 0
    )
    
    lambda <- slice_sampling(
      x0 = lambda,
      g = function(lambda) dgamma(
        x = lambda, shape = a_lambda, rate = b_lambda, log = TRUE
      ) + sum(dnbinom(x = N, size = r, mu = lambda, log = TRUE)),
      w = sqrt(a_lambda) / b_lambda,
      lower = 0
    )
    
    
    if (iter > burnin & (iter - burnin) %% thin == 0L) {
      
      it <- (iter - burnin) / thin
      
      beta_Gibbs[it, ] <- beta
      gamma_Gibbs[it, ] <- gamma
      
      N_Gibbs[it, ] <- N
      Y_unobs_Gibbs[it, , ] <- Y_unobs[which(cens), ]
      S_Gibbs[it, ] <- S
      
      s_Gibbs[it, ] <- s
      K_Gibbs[it] <- K
      m_Gibbs[it, , ] <- m
      delta_Gibbs[it, ] <- delta
      
      M_Gibbs[it] <- M
      
      sigma2_Gibbs[it] <- sigma2
      eta2_Gibbs[it] <- eta2
      r_Gibbs[it] <- r
      lambda_Gibbs[it] <- lambda
      
    }
    
    
    setTxtProgressBar(pb, iter)
  }
  
  close(pb)
  
  
  list(
    beta = beta_Gibbs,
    gamma = gamma_Gibbs,
    N = N_Gibbs,
    Y_unobs = Y_unobs_Gibbs,
    S = S_Gibbs,
    s = s_Gibbs,
    K = K_Gibbs,
    m = m_Gibbs,
    delta = delta_Gibbs,
    M = M_Gibbs,
    sigma2 = sigma2_Gibbs,
    eta2 = eta2_Gibbs,
    r = r_Gibbs,
    lambda = lambda_Gibbs
  )
  
}