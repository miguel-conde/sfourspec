library(dplyr)

# https://www.coinbase.com/


spectral_analysis <- function(x) {

  N <- length(x)
  w_0 <- 2*pi/N

  x_k <- fft(x)

  P_avg <- mean(x^2)   # Potencia media = Valor cuadrático medio
  P_dc <- mean(x)      # Potencia media de la componente continua  - BIAS
  P_ac <- P_avg - P_dc # Potencia media de la componente alterna = varianza de x
  x_ef <- sqrt(P_ac)   # Valor eficaz

  tb_spec = tibble(n = 0:(N-1),
                   k = calc_w_k(N),
                   w_k = k*w_0,
                   f_k = w_k/2/pi,
                   T_k = 1/f_k,
                   x_k = x_k,
                   a_k = x_k/N,
                   # Th. Parseval: sum(F_L_spectrum) = P_avg = valor cuadrático medio
                   F_L_spectrum = Mod(a_k)^2)

  tb_dens_spec <- tb_spec %>%
    select(k, w_k, f_k, T_k, a_k, F_L_spectrum) %>%
    group_by(k, w_k, f_k, T_k) %>%
    summarise(periodogram_k = sum(F_L_spectrum),
              B_k = sum(a_k),
              C_k = ifelse(length(a_k) == 2, (0+1i)*(a_k[1]-a_k[2]), a_k)) %>%
    ungroup() %>%
    mutate(A_k = sqrt(Mod(B_k)^2 + Mod(C_k)^2))

  out <- list(P_avg = P_avg,
              P_dc = P_dc,
              P_ac = P_ac,
              x_ef = x_ef,
              tb_spec = tb_spec,
              tb_dens_spec = tb_dens_spec)
  return(out)
}

calc_w_k <- function(N) {

  out <- rep(NA, N)
  out[1] <- 0

  if (N %% 2 == 0) { # N par

    for (i in 1:(N/2-1)) {
      out[i+1] <- out[N-i+1] <- i
    }

    out[N/2+1] <- N/2

  } else { # N impar

    for (i in 1:(N/2)) {
      out[i+1] <- out[N-i+1] <- i
    }
  }


  return(out)
}

rebuild_signal <- function(spec_x, threshold = .8) {
  # Armónicos menos potentes
  aux <- spec_x$tb_dens_spec %>%
    arrange(-periodogram_k)
  ks <- aux %>%
    .$k %>% .[which(cumsum(aux$periodogram_k) > threshold*spec_x$P_ac)]

  # Reconstruimos la señal solo con los armónicos más potentes
  aux2 <-  spec_x$tb_spec %>%
    pull(a_k)
  aux2[ks] <- 0

  hat_x_season <- aux2 %>%
    fft(inverse = TRUE) %>%
    Re()

  return(hat_x_season)
}

rebuild_signal_harmonics <- function(spec_x, Ts) {
  # Armónicos deseados
  ks <- spec_x$tb_dens_spec %>%
    filter(T_k %in% Ts) %>%
    pull(k)

  # Reconstruimos la señal solo con los armónicos deseados
  aux <-  spec_x$tb_spec %>%
    pull(a_k)
  aux[!(spec_x$tb_spec$k %in% ks)] <- 0

  hat_x_season <- aux %>%
    fft(inverse = TRUE) %>%
    Re()

  return(hat_x_season)
}