# depends: 
sir_piv_pit_feasible <- function(PIV, PIT, s0, i0, rtol = 1e-8) {
  # Test whether a proposed (PIV, PIT) pair is achievable by a standard
  # continuous-time SIR model
  #
  #   ds/dt = -beta * s * i
  #   di/dt =  beta * s * i - gamma * i
  #
  # with instantaneous incidence
  #
  #   incidence(t) = beta * s(t) * i(t),
  #
  # for some beta > 0 and gamma > 0, with fixed initial conditions s0, i0.
  #
  # Parameters
  # ----------
  # PIV  : proposed peak incidence value
  # PIT  : proposed peak incidence time
  # s0   : initial susceptible population proportion
  # i0   : initial infected population proportion
  # rtol : numerical relative tolerance
  #
  # Returns
  # -------
  # TRUE if the pair is numerically feasible, FALSE otherwise.

  # Basic checks
  if (PIV < 0 || PIT < 0) {
    return(FALSE)
  }

  if (s0 <= 0 || i0 <= 0) {
    return(FALSE)
  }

  if (s0 + i0 > 1 + rtol) {
    return(FALSE)
  }

  # PIT = 0:
  #
  # If the incidence maximum occurs initially, any positive PIV can
  # be obtained by scaling beta, provided incidence can be made
  # non-increasing initially.
  if (PIT == 0) {
    return(PIV > 0)
  }

  # A positive-time incidence peak requires s0 - i0 > 0.
  kcrit <- s0 - i0

  if (kcrit <= 0) {
    return(FALSE)
  }

  if (PIV <= 0) {
    return(FALSE)
  }

  target <- PIV * PIT

  # H(k) = dimensionless peak incidence * dimensionless peak time
  #
  # k = gamma / beta
  # Positive-time peaks require 0 < k < s0 - i0.
  H <- function(k) {
    q <- s0 + i0

    # At the incidence peak:
    #
    #   s* - i* = k
    #
    # together with the SIR invariant gives
    #
    #   2 s* - q - k [1 + log(s*/s0)] = 0.
    peak_equation <- function(s) {
      2 * s - q - k * (1 + log(s / s0))
    }

    root <- uniroot(
      peak_equation,
      interval = c(k, s0),
      tol = 1e-12
    )

    s_star <- root$root
    i_star <- s_star - k
    m <- s_star * i_star

    # Along the SIR trajectory:
    #
    #   i(s) = s0 + i0 - s + k log(s/s0)
    i_of_s <- function(s) {
      q - s + k * log(s / s0)
    }

    # Dimensionless time to incidence peak, u = beta * t
    integrand <- function(s) {
      1 / (s * i_of_s(s))
    }

    B <- integrate(
      integrand,
      lower = s_star,
      upper = s0,
      rel.tol = 1e-10,
      abs.tol = 1e-10,
      subdivisions = 200L
    )$value

    m * B
  }

  # Include the k -> 0 limit analytically:
  #
  #   H(0) = (s0 + i0)/4 * log(s0/i0)
  H0 <- (s0 + i0) / 4 * log(s0 / i0)

  # Search for a possible interior maximum, staying slightly away
  # from the boundaries for numerical stability.
  eps <- 1e-9 * kcrit

  opt <- optimize(
    f = H,
    interval = c(eps, kcrit - eps),
    maximum = TRUE,
    tol = 1e-10
  )

  H_interior <- opt$objective
  Hmax <- max(H0, H_interior)

  # H(k) approaches 0 as k -> kcrit.
  #
  # If Hmax occurs only at k = 0, it corresponds to gamma = 0,
  # which is excluded here. Therefore equality at that upper
  # attainable-product boundary is treated as infeasible unless an
  # interior k attains the same value numerically.
  strictly_below <- target < Hmax

  approximately_equal <- abs(target - Hmax) <= rtol * abs(Hmax)

  interior_attains_max <- H_interior > H0 * (1 + rtol)

  strictly_below || (approximately_equal && interior_attains_max)
}

