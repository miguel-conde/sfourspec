# https://math.stackexchange.com/questions/239134/fourier-series-and-exponential
# https://lpsa.swarthmore.edu/Fourier/Series/DerFS.html#Equiv

#' spectral_analysis
#'
#' Simple Fourier spectra
#'
#' @param x a numeric vector or univariate time series.
#'
#' @return List with:
#'     \itemize{
#'         \item P_avg \code{double}, Mean Square Value of x, aka x's mean power.
#'         \item P_dc \code{double}, Square mean value of x, aka the mean power of the
#'         x's continuous component (dc).
#'         \item P_ac \code{double}, Variance of x, aka the mean power of the x's
#'         alternate component (ac).
#'         \item x_ef \code{double},, The effective value of x.
#'         \item four_exp \code{tibble}, containing the results obtained using a
#'         Fourier series with exponential form (linear combination of complex
#'         exponentials).
#'         \item four_cos_sin \code{tibble}, containing the results obtained using a
#'         Fourier series with sine / cosine form (linear combination of sines and cosines)..
#'         \item four_cos \code{tibble}, containing the results obtained using a
#'         Fourier series with cosine form (linear combination of cosines).
#'         \item x, the original signal.
#'     }
#'
#' @details
#'
#' \deqn{
#' \langle x[n] \rangle_{N} = \frac{1}{N} \sum_{\langle N \rangle}x[n]
#' }
#'
#' @export
#'
#' @importFrom stats acf fft
#'
#' @encoding UTF-8
#'
#' @examples
#' N_0 <- 10       # PERIODO - u. de t. (p. ej., segundos)
#' w_0 <- 2*pi/N_0 # FRECUENCIA ANGULAR
#'
#' # (muestreo cada milisegundo)
#' sampling_freq <- 1e-3
#'
#' ns <- seq(0, N_0/sampling_freq-1, by = 1)
#' N <- length(ns)
#' x <- sin(w_0 * ns * sampling_freq)
#'
#' spec_x <- spectral_analysis(x)
#'
#' plot(spec_x)
#'
spectral_analysis <- function(x) {

  N <- length(x)
  w_0 <- 2*pi/N

  x_k <- fft(x)

  P_avg <- mean(x^2)   # Potencia media = Valor cuadrático medio
  P_dc <- mean(x)      # Potencia media de la componente continua  - BIAS
  P_ac <- P_avg - P_dc # Potencia media de la componente alterna = varianza de x
  x_ef <- sqrt(P_ac)   # Valor eficaz

  four_exp = tibble::tibble(n = 0:(N-1),
                            k = calc_w_k(N),
                            w_k = k*w_0,
                            f_k = w_k/2/pi,
                            T_k = 1/f_k,
                            x_k = x_k,
                            a_k = x_k/N,
                            # Th. Parseval: sum(F_L_spectrum) = P_avg = valor cuadrático medio
                            F_L_spectrum_k = Mod(a_k)^2,
                            periodogram_k = calc_periodogram(x, P_avg))

  four_cos_sin <- four_exp %>%
    select(k, w_k, f_k, T_k, a_k, F_L_spectrum_k, periodogram_k) %>%
    group_by(k, w_k, f_k, T_k) %>%
    summarise(B_k = sum(a_k),
              C_k = ifelse(length(a_k) == 2, (0+1i)*(a_k[1]-a_k[2]), a_k),
              F_L_spectrum_k = sum(F_L_spectrum_k),
              periodogram_k = sum(periodogram_k)) %>%
    ungroup()

  four_cos <- four_cos_sin %>%
    mutate(A_k = sqrt(Mod(B_k)^2 + Mod(C_k)^2),
           theta = atan(-C_k / A_k),
           F_L_spectrum_k = 1/2*Mod(A_k)^2) %>%
    select(-B_k, -C_k)

  out <- list(P_avg = P_avg,
              P_dc = P_dc,
              P_ac = P_ac,
              x_ef = x_ef,
              four_exp = four_exp,
              four_cos_sin = four_cos_sin,
              four_cos = four_cos,
              x = x)

  class(out) <- c("sfourspec", class(out))

  return(out)
}

#' calc_periodogram
#'
#' This function estimates the spectral density as the Fourier Transform of the
#' auto covariance function.
#'
#' @param x a numeric vector or univariate time series.
#' @param P_avg  \code{double}, Mean Square Value of x, aka x's mean power.
#'
#' @return the spectral density estimation.
#'
#' @encoding UTF-8
#' @keywords internal
#'
#' @examples
#'
calc_periodogram <- function(x, P_avg) {

  N <- length(x)

  acf_x <- acf(x, type = "covariance", lag.max = N, plot = FALSE)$acf %>%
    as.numeric()

  fft_acf <- Mod(fft(acf_x))

  area_fft_acf <- mean(fft_acf)*0.5

  fft_acf <- fft_acf/area_fft_acf*P_avg

  fft_acf
}

#' calc_w_k
#'
#' @param N TODO
#'
#' @return
#'
#' @encoding UTF-8
#' @keywords internal
#'
#' @examples
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

#' rebuild_signal
#'
#' @param spec_x an object of class \code{sfourspec}.
#' @param threshold the % of power in the original signal to keep in the rebuilt
#' signal
#'
#' @return a numeric vector with the rebuilt signal.
#' @export
#'
#' @encoding UTF-8
#'
#' @examples
rebuild_signal <- function(spec_x, threshold = .8) {
  # Armónicos menos potentes
  aux <- spec_x$four_cos_sin %>%
    arrange(-F_L_spectrum_k)
  ks <- aux %>%
    .$k %>% .[which(cumsum(aux$F_L_spectrum_k) > threshold*spec_x$P_avg)]

  # Reconstruimos la señal solo con los armónicos más potentes
  aux2 <-  spec_x$four_exp %>%
    pull(a_k)
  aux2[ks] <- 0

  hat_x <- aux2 %>%
    fft(inverse = TRUE) %>%
    Re()

  if (inherits(spec_x$x, "ts")) {
    hat_x <- make_ts(hat_x, spec_x)
  }

  return(hat_x)
}

#' make_ts
#'
#' @param hat_x
#' @param spec_x
#'
#' @return
#'
#' @importFrom stats ts start end frequency deltat
#'
#' @encoding UTF-8
#' @keywords internal
#'
#' @examples
make_ts <- function(hat_x, spec_x) {

  out <- ts(hat_x,
              start = start(spec_x$x),
              end = end(spec_x$x),
              frequency = frequency(spec_x$x),
              deltat = deltat(spec_x$x),
              class = class(spec_x$x),
              names = names(spec_x$x))

  return(out)
}

#' rebuild_signal_harmonics
#'
#' @param spec_x TODO
#' @param Ts TODO
#'
#' @return
#' @export
#'
#' @encoding UTF-8
#'
#' @examples
rebuild_signal_harmonics <- function(spec_x, Ts) {
  # Armónicos deseados
  ks <- spec_x$four_cos_sin %>%
    filter(T_k %in% Ts) %>%
    pull(k)

  # Reconstruimos la señal solo con los armónicos deseados
  aux <-  spec_x$four_exp %>%
    pull(a_k)
  aux[!(spec_x$four_exp$k %in% ks)] <- 0

  hat_x <- aux %>%
    fft(inverse = TRUE) %>%
    Re()

  if (inherits(spec_x$x, "ts")) {
    hat_x <- make_ts(hat_x, spec_x)
  }

  return(hat_x)
}

#' plot_spectrum
#'
#' @param x TODO
#' @param fft_type TODO
#' @param x_axis TODO
#' @param y_axis TODO
#' @param ... TODO
#'
#' @return
#'
#' @export
#'
#' @importFrom graphics plot
#'
#' @encoding UTF-8
#'
#' @examples
plot.sfourspec <- function(x,
                           fft_type = "exp",
                           x_axis = "f",
                           y_axis = "fourier_line_spectrum",
                           ...) {

  spec_x <- x

  fft_type <- switch(match.arg(fft_type, c("exp", "cos_sin", "cos")),
                     "exp" = "four_exp",
                     "cos_sin" = "four_cos_sin",
                     "cos" = "four_cos")
  x_axis = switch(match.arg(x_axis, c("f", "w", "T", "k", "n")),
                  "f" = "f_k",
                  "w" = "w_k",
                  "T" = "T_k",
                  "k" = "k",
                  "n" = "n")
  y_axis = switch(match.arg(y_axis, c("fourier_line_spectrum",
                                 "periodogram")),
                  "fourier_line_spectrum" = "F_L_spectrum_k",
                  "periodogram" = "periodogram_k")

  if (x_axis == "n" & fft_type != "four_exp") {
    stop(paste("Can't use 'n' with fft_type =", fft_type))
  }

  title = switch(y_axis,
                 "F_L_spectrum_k" = "Fourier Line Spectrum",
                 "periodogram_k" = "Periodogram (density spectrum)")
  type_plot = switch(y_axis,
                 "F_L_spectrum_k" = "h",
                 "periodogram_k" = "l")
  x_lab = switch(x_axis,
                  "f_k" = expression(paste(f[k], " (u.de ", t^-1, ")")),
                  "w_k" = expression(paste(omega[k], " (cycles /u.de t.")),
                  "T_k" = expression(paste(T[k], " (u.de t.)")),
                  "k" = "k",
                  "n" = "n")
  y_lab = switch(y_axis,
                 "F_L_spectrum_k" = expression(paste("Var(x) = ", P[ac])),
                 "periodogram_k" = expression(paste(P[ac], " / u.de ", t^-1)))

# browser()
  plot(spec_x[[fft_type]] %>% select(all_of(c(x_axis, y_axis))),
       type = type_plot,
       xlab = x_lab,
       ylab = y_lab,
       main = title,
       ...)
}
