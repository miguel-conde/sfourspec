context("spectral_analysis()")

library(sfourspec)

test_that("The number of Fourier coefficients equals the number of samples",
          {
            air_pass <- tibble::tibble(t = 1:length(AirPassengers),
                                       pass = as.numeric(AirPassengers))

            lm_air_pass <- lm(pass ~ t, air_pass)

            detrended_air_pass <- air_pass$pass - fitted(lm_air_pass)

            N <- length(detrended_air_pass)
            ns <- 0:(N-1)

            spec_x <- spectral_analysis(detrended_air_pass)

            # Even sample
            expect_equal(nrow(spec_x$four_exp), N)

            spec_x_2 <- spectral_analysis(detrended_air_pass)[-1]

            # Odd sample
            expect_equal(nrow(spec_x_2$four_exp), N-1)
          })
