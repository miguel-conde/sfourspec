context("spectral_analysis()")

library(sfourspec)

N_0 <- 10       # PERIODO - u. de t. (p. ej., segundos)
w_0 <- 2*pi/N_0 # FRECUENCIA ANGULAR

# (muestreo cada milisegundo)
sampling_freq <- 1e-3

ns <- seq(0, N_0/sampling_freq-1, by = 1) # (muestreo cada milisegundo)
N <- length(ns)
x <- sin(w_0*ns*sampling_freq)

spec_x <- spectral_analysis(x)       # Even sample
spec_x_2 <- spectral_analysis(x[-1]) # Odd sample

test_that("The powers add up OK",
          {
            expect_equal(spec_x$P_avg, spec_x$P_dc + spec_x$P_ac)

          })

test_that("The Power_ac equals the theorical variance of the sine",
          {
            expect_equal(abs(spec_x$P_ac - 0.5), .Machine$double.neg.eps)

          })

test_that("The number of Fourier coefficients equals the number of samples",
          {
            # Even sample
            expect_equal(nrow(spec_x$four_exp), N)

            # Odd sample
            expect_equal(nrow(spec_x_2$four_exp), N-1)
          })

test_that("The number of armonics is ok",
          {
            expect_equal(range(spec_x$four_exp$k),
                         as.integer(range(spec_x$four_exp$n) / 2 + 0:1))

            expect_equal(range(spec_x_2$four_exp$k),
                         range(spec_x_2$four_exp$n) / 2)
          })

test_that("The sum of each Fourier Line Spectra adds up to the Mean Square Value",
          {
            expect_equal(abs(sum(spec_x$four_exp$F_L_spectrum_k) - spec_x$P_avg),
                         .Machine$double.neg.eps)
            expect_equal(abs(sum(spec_x$four_cos_sin$F_L_spectrum_k) - spec_x$P_avg),
                         .Machine$double.neg.eps)
            expect_equal(abs(sum(spec_x$four_cos$F_L_spectrum_k) - spec_x$P_avg),
                         .Machine$double.neg.eps)
          })

test_that("The sum of each Periodogram adds up to the number of samples",
          {
            expect_equal(abs(sum(spec_x$four_exp$periodogram_k)), N)
            expect_equal(abs(sum(spec_x$four_cos_sin$periodogram_k)), N)
            expect_equal(abs(sum(spec_x$four_cos$periodogram_k)), N)
          })
