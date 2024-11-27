# Modelo 1 - Dados Misra1a
Misra1a <- data.frame(y = c(10.07, 14.73, 17.94, 23.93, 29.61,
                            35.18, 40.02, 44.82, 50.76, 55.05,
                            61.01, 66.4, 75.47, 81.78),
                      x = c(77.6, 114.9, 141.1, 190.8, 239.9,
                            289, 332.8, 378.4, 434.8, 477.3,
                            536.8, 593.1, 689.1, 760))
beta_init1 <- c(500, 0.0001)  # Primeiro valor inicial
beta_init2 <- c(250, 0.0005)  # Segundo valor inicial

# Função para calcular o gradiente
gradient <- function(beta, data) {
  y <- data$y
  x <- data$x
  g1 <- sum((-2 + 2 * exp(-beta[2] * x)) * (-beta[1] * (1 - exp(-beta[2] * x)) + y))
  g2 <- sum(-2 * beta[1] * x * (-beta[1] * (1 - exp(-beta[2] * x)) + y) * exp(-beta[2] * x))
  return(c(g1, g2))
}

# Função para calcular a matriz Hessiana
hessian <- function(beta, data) {
  y <- data$y
  x <- data$x
  h11 <- sum((-2 + 2 * exp(-beta[2] * x)) * (-1 + exp(-beta[2] * x)))
  h12 <- sum(-beta[1] * x * (-2 + 2 * exp(-beta[2] * x)) * exp(-beta[2] * x) -
               2 * x * (-beta[1] * (1 - exp(-beta[2] * x)) + y) * exp(-beta[2] * x))
  h21 <- h12
  h22 <- sum(2 * beta[1]^2 * x^2 * exp(-2 * beta[2] * x) +
               2 * beta[1] * x^2 * (-beta[1] * (1 - exp(-beta[2] * x)) + y) * exp(-beta[2] * x))
  return(matrix(c(h11, h12, h21, h22), nrow = 2, byrow = TRUE))
}

# Função do método de Newton-Raphson 
newton_raphson <- function(beta_init, data, epsilon = 1e-6, max_iter = 1000, gamma = 0.4) {
  beta <- beta_init
  for (i in 1:max_iter) {
    grad <- gradient(beta, data)  
    hess <- hessian(beta, data)  
    
    hess_regularized <- hess + diag(gamma, nrow(hess))
    
    
    delta_beta <- solve(hess_regularized) %*% grad  
    
    beta_new <- beta - delta_beta
    
    if (sqrt(sum((beta_new - beta)^2)) < epsilon) {
      return(beta_new)
    }
    
    beta <- beta_new
  }
  
  return(beta)
}


# Resultados do método de Newton-Raphson
result1 <- newton_raphson(beta_init1, Misra1a)
result1

result2 <- newton_raphson(beta_init2, Misra1a)
result2


# Função para calcular a função objetivo
objective_function <- function(beta, data) {
  y <- data$y
  x <- data$x
  residuals <- y - beta[1] * (1 - exp(-beta[2] * x))
  return(sum(residuals^2))
}

# Função para criar o gráfico de line search
line_search_plot <- function(beta, grad, hess, data, gamma_candidates) {
  direction <- -solve(hess) %*% grad
  objective_values <- numeric(length(gamma_candidates))
  
  # Avalia a função objetivo para cada gamma
  for (i in seq_along(gamma_candidates)) {
    gamma <- gamma_candidates[i]
    beta_new <- beta + gamma * direction
    objective_values[i] <- objective_function(beta_new, data)
  }
  
  plot(gamma_candidates, objective_values, type = "b", col = "blue",
       main = "Line Search", xlab = expression(gamma),
       ylab = "Objective Function")
}

gamma_candidates <- seq(0.01, 1, by = 0.01)

# Gera o gráfico
line_search_plot(beta_init1, gradient(beta_init1, Misra1a), hessian(beta_init1, Misra1a), Misra1a, gamma_candidates)
line_search_plot(beta_init2, gradient(beta_init2, Misra1a), hessian(beta_init2, Misra1a), Misra1a, gamma_candidates)


# Função de line search
line_search <- function(beta, grad, hess, data, gamma_candidates) {
  direction <- -solve(hess) %*% grad
  for (gamma in gamma_candidates) {
    beta_new <- beta + gamma * direction
    if (objective_function(beta_new, data) < objective_function(beta, data)) {
      return(list(beta_new = beta_new, gamma = gamma))
    }
  }
  return(list(beta_new = beta, gamma = NA))
}

# Método de Newton-Raphson com line search ajustado
newton_raphson_ls <- function(beta_init, data, epsilon = 1e-6, max_iter = 10000, gamma_candidates) {
  beta <- beta_init
  for (i in 1:max_iter) {
    grad <- gradient(beta, data)
    hess <- hessian(beta, data)
    
    # Verifica invertibilidade da Hessiana
    if (det(hess) == 0) {
      stop("A matriz Hessiana não é invertível.")
    }
    
    # Line search
    ls_result <- line_search(beta, grad, hess, data, gamma_candidates)
    beta_new <- ls_result$beta_new
    gamma <- ls_result$gamma
    
    cat("Iteração", i, "- Gamma:", gamma, "- Beta:", beta_new, "\n")
    
    # Critério de parada
    if (sqrt(sum((beta_new - beta)^2)) < epsilon) {
      cat("Convergência atingida na iteração", i, "\n")
      return(list(beta = beta_new, converged = TRUE))
    }
    
    beta <- beta_new
  }
  cat("Máximo de iterações alcançado sem convergência.\n")
  return(list(beta = beta, converged = FALSE))
}


# Resultados para o primeiro conjunto de valores iniciais
result1 <- newton_raphson_ls(beta_init1, Misra1a, gamma_candidates = 0.1)
result1

# Resultados para o segundo conjunto de valores iniciais
result2 <- newton_raphson_ls(beta_init2, Misra1a, gamma_candidates = 0.5)
result2



# Função de descida de gradiente 
gradient_descent_ls <- function(beta_init, data, gamma = 0.1, epsilon = 1e-6, max_iter = 1000) {
  beta <- beta_init
  for (k in 1:max_iter) {
    grad <- gradient(beta, data)
    gamma_k <- gamma * 2^(-k)
    beta_new <- beta - gamma_k * grad
    
    # Exibe a iteração atual
    cat("Iteração", k, "- Gamma_k:", gamma_k, "- Beta:", beta_new, "\n")
    
    # Critério de parada
    if (sqrt(sum((beta_new - beta)^2)) < epsilon) {
      cat("Convergência atingida na iteração", k, "\n")
      return(beta_new)
    }
    
    beta <- beta_new
  }
  cat("Máximo de iterações alcançado sem convergência.\n")
  return(beta)
}

result1 <- gradient_descent_ls(beta_init1, Misra1a)
result1

result2 <- gradient_descent_ls(beta_init2, Misra1a)
result2


# Função de descida de gradiente com gamma fixo
gradient_descent_fixed <- function(beta_init, data, gamma = 0.02, epsilon = 1e-6, max_iter = 1000) {
  beta <- beta_init
  for (k in 1:max_iter) {
    grad <- gradient(beta, data)
    beta_new <- beta - gamma * grad
    beta_new <- pmax(beta_new, 0)
    
    # Exibe a iteração atual
    cat("Iteração", k, "- Gamma:", gamma, "- Beta:", beta_new, "\n")
    
    # Critério de parada
    if (sqrt(sum((beta_new - beta)^2)) < epsilon) {
      cat("Convergência atingida na iteração", k, "\n")
      return(beta_new)
    }
    beta <- beta_new
  }
  
  cat("Máximo de iterações alcançado sem convergência.\n")
  return(beta)
}

result1 <- gradient_descent_fixed(beta_init1, Misra1a)
result1

result2 <- gradient_descent_fixed(beta_init2, Misra1a)
result2



#Modelo 2 - Dados Thurber
Thurber <- data.frame(y = c(80.574, 84.248, 87.264, 87.195, 89.076, 89.608,
                            89.868, 90.101, 92.405, 95.854, 100.696, 101.06,
                            401.672, 390.724, 567.534, 635.316, 733.054,
                            759.087, 894.206, 990.785, 1090.109, 1080.914,
                            1122.643, 1178.351, 1260.531, 1273.514, 1288.339,
                            1327.543, 1353.863, 1414.509, 1425.208, 1421.384,
                            1442.962, 1464.35, 1468.705, 1447.894, 1457.628),
                      x = c(-3.067, -2.981, -2.921, -2.912, -2.84, -2.797,
                            -2.702, -2.699, -2.633, -2.481, -2.363, -2.322,
                            -1.501, -1.46, -1.274, -1.212, -1.1, -1.046,
                            -0.915, -0.714, -0.566, -0.545, -0.4, -0.309,
                            -0.109, -0.103, 0.01, 0.119, 0.377, 0.79, 0.963,
                            1.006, 1.115, 1.572, 1.841, 2.047, 2.2))

beta0_chapeu_thurber_1 <- c(1000, 1000, 400, 40, 0.7, 0.3, 0.03)
beta0_chapeu_thurber_2 <- c(1300, 1500, 500, 75, 1, 0.4, 0.05)

# Gradiente
gradient_thurber <- function(beta, data) {
  x <- data$x
  y <- data$y
  
  num <- beta[1] + beta[2] * x + beta[3] * x^2 + beta[4] * x^3
  denom <- 1 + beta[5] * x + beta[6] * x^2 + beta[7] * x^3
  f <- num / denom
  
  grad <- numeric(7)
  grad[1] <- sum(-2 * (y - f) / denom)
  grad[2] <- sum(-2 * x * (y - f) / denom)
  grad[3] <- sum(-2 * x^2 * (y - f) / denom)
  grad[4] <- sum(-2 * x^3 * (y - f) / denom)
  grad[5] <- sum(2 * x * (y - f) * (num / denom^2))
  grad[6] <- sum(2 * x^2 * (y - f) * (num / denom^2))
  grad[7] <- sum(2 * x^3 * (y - f) * (num / denom^2))
  
  return(grad)
}

#Hessiana
hessian_thurber <- function(beta, data) {
  x <- data$x
  y <- data$y
  num <- beta[1] + beta[2] * x + beta[3] * x^2 + beta[4] * x^3
  denom <- 1 + beta[5] * x + beta[6] * x^2 + beta[7] * x^3
  f <- num / denom
  hessian <- matrix(0, nrow = 7, ncol = 7)
  
  for (i in 1:length(x)) {

    common_term <- y[i] - f[i]
    d1 <- 1 / denom[i]
    d2 <- num[i] / denom[i]^2
    d3 <- denom[i]^3
    d4 <- denom[i]^4
    
    # Diagonais principais
    hessian[1, 1] <- hessian[1, 1] + 2 * d1^2
    hessian[2, 2] <- hessian[2, 2] + 2 * (x[i] * d1)^2
    hessian[3, 3] <- hessian[3, 3] + 2 * (x[i]^2 * d1)^2
    hessian[4, 4] <- hessian[4, 4] + 2 * (x[i]^3 * d1)^2
    hessian[5, 5] <- hessian[5, 5] + 2 * x[i]^2 * ((-4 * common_term * d2) + (num[i] * x[i]^2 / d4))
    hessian[6, 6] <- hessian[6, 6] + 2 * x[i]^4 * ((-4 * common_term * d2) + (num[i] * x[i]^4 / d4))
    hessian[7, 7] <- hessian[7, 7] + 2 * x[i]^6 * ((-4 * common_term * d2) + (num[i] * x[i]^6 / d4))
    
    # Termos cruzados
    for (j in 1:7) {
      for (k in j:7) {
        if (j <= 4 && k <= 4) {
          hessian[j, k] <- hessian[j, k] + 2 * (x[i]^(j - 1) * x[i]^(k - 1) * d1^2)
        } else if (j > 4 || k > 4) {
          # Interações cruzadas com beta_5, beta_6, beta_7
          hessian[j, k] <- hessian[j, k] + 
            2 * (x[i]^(j - 4) * x[i]^(k - 4) * 
                   ((common_term * d2) - (num[i] * x[i]^(j + k - 8) / d4)))
        }
      }
    }
  }
  
  hessian[lower.tri(hessian)] <- t(hessian)[lower.tri(hessian)]
  
  return(hessian)
}


# Função do método de Newton-Raphson 
newton_raphson_thurber <- function(beta_init, data, epsilon = 1e-6, max_iter = 10000, gamma = 1, max_gamma_iter = 10) {
  x <- data$x
  y <- data$y
  
  beta <- beta_init  
  
  for (i in 1:max_iter) {
    grad <- gradient_thurber(beta, data)
    hess <- hessian_thurber(beta, data)
    
    gamma_adjust <- gamma
    success <- FALSE
    for (gamma_iter in 1:max_gamma_iter) {
      hess_regularized <- hess + diag(gamma_adjust, nrow(hess))
      eigenvalues <- eigen(hess_regularized)$values
      if (all(is.finite(eigenvalues)) && min(abs(eigenvalues)) > 1e-10) {
        success <- TRUE
        break  # Hessiana regularizada é aceitável
      }
      gamma_adjust <- gamma_adjust * 10  # Aumenta gamma
    }
    
    if (!success) {
      cat("Hessiana regularizada ainda singular na iteração", i, "\n")
      return(beta)
    }
    
    # Calcula o passo de atualização
    delta_beta <- tryCatch({
      solve(hess_regularized, grad)
    }, error = function(e) {
      cat("Erro ao resolver o sistema na iteração", i, "\n")
      return(rep(NA, length(beta)))
    })
    
    # Verifica se delta_beta é válido
    if (any(is.na(delta_beta))) {
      cat("Delta beta inválido na iteração", i, "\n")
      return(beta) 
    }
    
    # Reduz o tamanho do passo se necessário
    step_size <- sqrt(sum(delta_beta^2))
    if (step_size > 10) {  
      delta_beta <- delta_beta / step_size
    }
    
    beta_new <- beta - delta_beta
    
    # Critério de parada
    if (sqrt(sum((beta_new - beta)^2)) < epsilon) {
      cat("Convergência atingida na iteração", i, "\n")
      return(beta_new)
    }
    
    beta <- beta_new
  }
  
  cat("Máximo de iterações alcançado sem convergência.\n")
  return(beta)
}

resultado1 <- newton_raphson_thurber(beta0_chapeu_thurber_1, Thurber)
resultado1

resultado2 <- newton_raphson_thurber(beta0_chapeu_thurber_2, Thurber)
resultado2


# Função objetivo
objective_function <- function(beta, data) {
  y <- data$y
  x <- data$x
  num <- beta[1] + beta[2] * x + beta[3] * x^2 + beta[4] * x^3
  denom <- 1 + beta[5] * x + beta[6] * x^2 + beta[7] * x^3
  residuals <- y - (num / denom)
  return(sum(residuals^2))
}

# Função para realizar line-search utilizando um grid de valores para γ
line_search <- function(beta, grad, hess, data, gamma_candidates) {
  direction <- -solve(hess) %*% grad
  best_gamma <- NA
  best_objective <- Inf
  
  for (gamma in gamma_candidates) {
    beta_new <- beta + gamma * direction
    obj_value <- objective_function(beta_new, data)
    if (obj_value < best_objective) {
      best_objective <- obj_value
      best_gamma <- gamma
    }
  }
  
  return(best_gamma)
}

# Função de Newton-Raphson com Line Search 
newton_raphson_line_search <- function(beta_init, data, epsilon = 1e-6, max_iter = 10000, gamma_candidates = seq(0.01, 0.1, by = 0.01)) {
  beta <- beta_init
  objective_values <- c()  
  beta_history <- matrix(NA, nrow = max_iter, ncol = length(beta_init))  
  
  for (i in 1:max_iter) {
    grad <- gradient_thurber(beta, data)
    hess <- hessian_thurber(beta, data)
    
    # Regularização robusta da Hessiana
    hess_regularized <- hess + diag(max(abs(hess)) * 0.1, nrow(hess))
    
    # Verifica se a Hessiana regularizada é invertível
    if (det(hess_regularized) == 0) {
      cat("Hessiana singular na iteração", i, "\n")
      break
    }
    
    # Line-search adaptativo
    gamma <- line_search(beta, grad, hess_regularized, data, gamma_candidates)
    direction <- -solve(hess_regularized) %*% grad
    beta_new <- beta + gamma * direction
    
    # Limita o tamanho do passo
    step_size <- sqrt(sum(direction^2))
    if (step_size > 1) {
      beta_new <- beta + (direction / step_size)
    }
    obj_value <- objective_function(beta_new, data)
    objective_values <- c(objective_values, obj_value)
    beta_history[i, ] <- beta_new
    
    # Critério de parada
    if (sqrt(sum((beta_new - beta)^2)) < epsilon) {
      cat("Convergência atingida na iteração", i, "\n")
      return(list(beta = beta_new, objective_values = objective_values, beta_history = beta_history[1:i, ]))
    }
    
    beta <- beta_new
  }
  
  cat("Máximo de iterações alcançado sem convergência.\n")
  return(list(beta = beta, objective_values = objective_values, beta_history = beta_history))
}


resultado1 <- newton_raphson_line_search(beta0_chapeu_thurber_1, Thurber)
resultado1$beta

resultado2 <- newton_raphson_line_search(beta0_chapeu_thurber_2, Thurber)
resultado2$beta



#Gráfico de convergência
plot(resultado1$objective_values, type = "l", col = "blue", xlab = "Iteração",
     ylab = "Função Objetivo", main = "Convergência - Beta Init 1")
lines(resultado2$objective_values, col = "red")
legend("topright", legend = c("Beta Init 1", "Beta Init 2"), col = c("blue", "red"), lty = 1)



# Função para realizar line-search adaptativo
line_search <- function(beta, grad, data, direction, gamma_candidates) {
  best_gamma <- NA
  best_objective <- Inf
  
  for (gamma in gamma_candidates) {
    beta_new <- beta + gamma * direction
    obj_value <- objective_function(beta_new, data)
    if (obj_value < best_objective) {
      best_objective <- obj_value
      best_gamma <- gamma
    }
  }
  
  return(best_gamma)
}


# Função de Normalização dos Dados
normalize_data <- function(data) {
  data$x <- scale(data$x)
  data$y <- scale(data$y)
  return(data)
}

# Função Objetivo Regularizada
objective_function_regularized <- function(beta, data, lambda = 1e-3) {
  y <- data$y
  x <- data$x
  num <- beta[1] + beta[2] * x + beta[3] * x^2 + beta[4] * x^3
  denom <- 1 + beta[5] * x + beta[6] * x^2 + beta[7] * x^3
  residuals <- y - (num / denom)
  return(sum(residuals^2) + lambda * sum(beta^2))  # Penalização L2
}

# Busca de Linha com Grid
line_search <- function(beta, grad, data, direction, gamma_candidates, lambda = 1e-3) {
  best_gamma <- NA
  best_objective <- Inf
  for (gamma in gamma_candidates) {
    beta_new <- beta + gamma * direction
    obj_value <- objective_function_regularized(beta_new, data, lambda)
    if (obj_value < best_objective) {
      best_objective <- obj_value
      best_gamma <- gamma
    }
  }
  return(best_gamma)
}

# # Função de Gradient Descent com Line-Search adaptativo
gradient_descent_line_search <- function(beta_init, data, gamma_candidates = seq(0.001, 0.2, by = 0.001),
                                         epsilon = 1e-6, max_iter = 10000, lambda = 1e-3) {
  beta <- beta_init
  objective_values <- c()
  beta_history <- matrix(NA, nrow = max_iter, ncol = length(beta_init))
  
  for (i in 1:max_iter) {
    grad <- gradient_thurber(beta, data)
    direction <- -grad  
    
    gamma <- line_search(beta, grad, data, direction, gamma_candidates, lambda)
    if (is.na(gamma)) {
      cat("Não foi possível encontrar um tamanho de passo adequado na iteração", i, "\n")
      break
    }
    
    # Atualização dos Parâmetros
    beta_new <- beta + gamma * direction
    obj_value <- objective_function_regularized(beta_new, data, lambda)
    objective_values <- c(objective_values, obj_value)
    beta_history[i, ] <- beta_new
    
    # Monitoramento do Gradiente
    cat(sprintf("Iteração %d: Norma do Gradiente = %.5f\n", i, sqrt(sum(grad^2))))
    
    # Critério de Parada
    if (sqrt(sum((beta_new - beta)^2)) < epsilon) {
      cat("Convergência atingida na iteração", i, "\n")
      return(list(beta = beta_new, objective_values = objective_values, beta_history = beta_history[1:i, ]))
    }
    
    beta <- beta_new
  }
  
  cat("Máximo de iterações alcançado sem convergência.\n")
  return(list(beta = beta, objective_values = objective_values, beta_history = beta_history))
}

# Normalização dos Dados
Thurber_normalized <- normalize_data(Thurber)

resultado1 <- gradient_descent_line_search(beta0_chapeu_thurber_1, Thurber_normalized)
resultado2 <- gradient_descent_line_search(beta0_chapeu_thurber_2, Thurber_normalized)

# Gráfico de Convergência
plot(resultado1$objective_values, type = "l", col = "blue", xlab = "Iteração",
     ylab = "Função Objetivo Regularizada", main = "Convergência - Beta Init 1 e Beta Init 2")
lines(resultado2$objective_values, col = "red")
legend("topright", legend = c("Beta Init 1", "Beta Init 2"), col = c("blue", "red"), lty = 1)


#Gráfico do movimento da solução
beta_history1 <- resultado1$beta_history
beta_history2 <- resultado2$beta_history

plot(Thurber$x, Thurber$y, pch = 20, col = "gray", main = "Movimento da Solução (Modelo 2)",
     xlab = "β1", ylab = "β2")
lines(beta_history1[, 1], beta_history1[, 2], type = "b", col = "blue", lty = 1)
points(beta_history2[, 1], beta_history2[, 2], type = "b", col = "red", pch = 19)
legend("topright", legend = c("Beta Init 1", "Beta Init 2"), col = c("blue", "red"), lty = 1)




