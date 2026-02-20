module modelutils

  use datadec

  interface
     function func(nparam, param, inc)
       real (kind=8) :: func
       integer (kind=4), intent(in) :: nparam
       real (kind=8), intent(in) :: param(*), inc
     end function func
  end interface

  type func_holder
     procedure(func), pointer, nopass :: f_ptr => null()
  end type func_holder

contains
  subroutine incdism (seed, nparam, param, incmin, incmax, inc, &
       dist, ierr, func)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine draws randomly a variable according to probability
! density \verb|func| with parameters \verb|param|. Same as previous
! routine, but can remember up to 10 different distributions at a time.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : October 2006
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     nparam: Number of parameters (I4)
!     param : Parameters (n*R8)
!     incmin: Minimum inclination (R8)
!     incmax: Maximum inclination (R8)
!     dist  : index of the selected distribution (I4)
!     func  : probability density function
!
! OUTPUT
!     inc   : Inclination (R8)
!     ierr  : Error code
!                0 : nominal run
!               10 : wrong input data
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) seed
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: param
!f2py intent(in) incmin
!f2py intent(in) incmax
!f2py intent(in) dist
!f2py intent(out) inc
!f2py intent(out) ierr
    implicit none

    integer (kind=4), intent(in) :: nparam, dist
    integer (kind=4), intent(inout) :: seed
    integer (kind=4), intent(out) :: ierr
    real (kind=8), intent(in) :: param(nparmax), incmin, incmax
    real (kind=8), intent(out) :: inc
    type(func_holder), intent(in) :: func
    integer :: i, ilo, ihi, di
    real (kind=8) :: random
    integer (kind=4), parameter :: np = 16384, nd = 10
    real (kind=8), save :: proba(0:np,nd), inctab(0:np,nd)
    logical, save :: first(nd)

    data first /.true.,.true.,.true.,.true.,.true., &
         .true.,.true.,.true.,.true.,.true./

    ierr = 0
    di = min(nd, dist)
    if (first(di)) then
       inctab(0,di) = incmin
       proba(0,di) = 0.d0
       do i = 1, np
          inctab(i,di) = incmin + dfloat(i)*(incmax-incmin)/dfloat(np)
          proba(i,di) = func%f_ptr(nparam, param(1:nparam), inctab(i,di)) + proba(i-1,di)
       end do
       do i = 1, np
          proba(i,di) = proba(i,di)/proba(np,di)
       end do
       first(di) = .false.
    end if

    random = ran_3(seed)
    inc = interp(proba(0,di), inctab(0,di), random, np+1)

    return
  end subroutine incdism

  real (kind=8) function Variably_tapered(h, params)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the number of objects brighter or equal to H
! following an exponentially tapered exponential with parameters in
! params.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : October 2021
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     params: parameters for the distribution (4*R8)
!
! OUTPUT
!     Variably_tapered : Random value of H (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) h
!f2py intent(in) params
    implicit none

    real (kind=8), intent(in) :: params(4), h

    Variably_tapered = 10.d0**(params(3)*3.d0*(h-params(1))/5.d0) &
         *exp(-10.d0**(-params(4)*3.d0*(h-params(2))/5.d0))

    return
  end function Variably_tapered

  real (kind=8) function Variably_tapered_diff(h, params)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the differential number of objects at mag H per unit
! mag following an exponentially tapered exponential with parameters in
! params.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : August 2023
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     params: parameters for the distribution (4*R8)
!
! OUTPUT
!     Variably_tapered_diff : Random value of H (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) h
!f2py intent(in) params
    implicit none

    real (kind=8), intent(in) :: params(4), h

    Variably_tapered_diff = (log(10.0d0)*3.0d0/5.0d0) &
         *10.0d0**(params(3)*3.0d0*(h-params(1))/5.0d0) &
         *(params(3) + params(4)*10.0d0**(-params(4)*3.0d0*(h-params(2))/5.0d0)) &
         *exp(-10.0d0**(-params(4)*3.0d0*(h-params(2))/5.0d0))

    return
  end function Variably_tapered_diff

  real (kind=8) function H_dist_cold(seed, h_max)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine draws randomly a number according to the cold belt H_r
! distribution, represented by an exponentially tapered exponential,
! with parameters I've fitted on the OSSOS cold belt data.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : October 2021
! Version 2 : December 2021- updated parameter values.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     h_max : Maximum value of H (R8)
!
! OUTPUT
!     H_dist_cold : Random value of H (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) seed
    implicit none

    integer (kind=4), parameter :: np = 16384
    integer (kind=4), intent(inout) :: seed
    real (kind=8), intent(in) :: h_max
    integer (kind=4) :: nparam, i
    real (kind=8) :: params(4), random, h_min
    real (kind=8), save :: proba(0:np), htab(0:np)
    logical, save :: first

    data first /.true./
    data params /-2.6d0, 8.1d0, 0.666d0, 0.42d0/
!    data params /-2.466d0, 7.895d0, 0.667d0, 0.438d0/
    data h_min /4.6d0/

    if (first) then
       htab(0) = h_min
       proba(0) = 1.d-10
       do i = 1, np
          htab(i) = h_min + dfloat(i)*(h_max-h_min)/dfloat(np)
          proba(i) = Variably_tapered(htab(i), params)
       end do
       do i = 0, np
          proba(i) = proba(i)/proba(np)
       end do
       first = .false.
    end if

    random = ran_3(seed)
    H_dist_cold = interp(proba, htab, random, np+1)

    return
  end function H_dist_cold

  real (kind=8) function H_dist_cold_2(seed, nparam, hparam)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine draws randomly a number according to the cold belt H_r
! distribution, represented by an exponentially tapered exponential,
! with parameters I've fitted on the OSSOS cold belt data.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : October 2021
! Version 2 : December 2021- updated parameter values.
! Version 3 : January 2022 - forcing slope at small sizes.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     nparam: Number of parameters (I4)
!     hparam: Parameters for asymptotic slope(s) (n*R8)
!             hparam(1): start of asymptote
!             hparam(2): contrast at start of asymptote
!             hparam(3): slope of asymptote
!             hparam(4): end of asymptote
!
! OUTPUT
!     H_dist_cold : Random value of H (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) seed
    implicit none

    integer (kind=4), parameter :: np = 16384
    integer (kind=4), intent(inout) :: seed
    integer (kind=4), intent(in) :: nparam
    real (kind=8), intent(in) :: hparam(4)
    integer (kind=4) :: i
    real (kind=8) :: params(4), random, h_min, h_max, n, c
    real (kind=8), save :: proba(0:np), htab(0:np)
    logical, save :: first

    data first /.true./
    data params /-2.6d0, 8.1d0, 0.666d0, 0.42d0/
    data h_min /4.6d0/

    if (first) then
       h_max = hparam(nparam)
       htab(0) = h_min
       proba(0) = 1.d-10
       n = Variably_tapered(hparam(1), params)
       c = hparam(2)
       do i = 1, np
          htab(i) = h_min + dfloat(i)*(h_max-h_min)/dfloat(np)
          if (htab(i) .lt. hparam(1)) then
             proba(i) = Variably_tapered(htab(i), params)
          else
             proba(i) = n &
                  + n*c*(10.0d0**(hparam(3)*(htab(i)-hparam(1))) - 1.0d0)
            end if
       end do
       do i = 0, np
          proba(i) = proba(i)/proba(np)
       end do
       first = .false.
    end if

    random = ran_3(seed)
    H_dist_cold_2 = interp(proba, htab, random, np+1)

    return
  end function H_dist_cold_2

  real (kind=8) function H_dist_hot(seed, h_max)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine draws randomly a number according to the hot belt H_r
! distribution, represented by an exponentially tapered exponential,
! with parameters I've fitted on the OSSOS hot belt data, in range 6-8.5.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : October 2021
! Version 2 : December 2021- updated parameter values.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     h_max : Maximum value of H (R8)
!
! OUTPUT
!     H_dist_hot : Random value of H (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) seed
    implicit none

    integer (kind=4), parameter :: np = 16384
    integer (kind=4), intent(inout) :: seed
    real (kind=8), intent(in) :: h_max
    integer (kind=4) :: nparam, i
    real (kind=8) :: params(4), random, h_min
    real (kind=8), save :: proba(0:np), htab(0:np)
    real (kind=8) :: n1, n2, n3, sl1, sl2, sl3, c1, cb, c, h1, h2, h3
    logical, save :: first

    data first /.true./
    data params /-2.465d0, 7.114d0, 0.666d0, 0.785d0/
    data h_min /-1.d0/
    data sl1 /0.13d0/, sl2 /0.6d0/, sl3 /0.4d0/, c /1.d0/
    data h1 /3.2d0/, h2 /6.d0/, h3 /8.5d0/, n1 /3.d0/

    if (first) then
       htab(0) = h_min
       proba(0) = 1.d-10
       n2 = Variably_tapered(h2, params)
       n3 = Variably_tapered(h3, params)
! The normalisation is done with the exponentially tapered exponential,
! as fitted on the OSSOS hot component. This determines the normalisation
! of the exponential between H = h2 and H = h1
! N(<H) = n2*10**(sl2*(H-h2))
! Then, there is an excess divot at h1. See
! [[file:///home/petit/Research/OSSOS/tes/OSSOSpapers/Papers/GlobalLuminosityFunction/CumDiffDistributions.org]]
! for the appropriate formula.
       cb = n1*sl1*log(10.d0)
       c1 = (n2-n1)*sl2/(n1*sl1*(10.d0**(sl2*(h2-h1))-1.d0))
       do i = 1, np
          htab(i) = h_min + dfloat(i)*(h_max-h_min)/dfloat(np)
          if (htab(i) .lt. h1) then
             proba(i) = n1*10.d0**(sl1*(htab(i)-h1))
          else if (htab(i) .lt. h2) then
             proba(i) = n1 + c1*cb &
                  *(10.d0**(sl2*(htab(i)-h1))-1.d0)/(sl2*log(10.d0))
          else if (htab(i) .lt. h3) then
             proba(i) = Variably_tapered(htab(i), params)
          else
             proba(i) = n3 + n3*c*(10.d0**(sl3*(htab(i)-h3)) - 1.d0)
          end if
       end do
       do i = 0, np
          proba(i) = proba(i)/proba(np)
       end do
       first = .false.
    end if

    random = ran_3(seed)
    H_dist_hot = interp(proba, htab, random, np+1)

    return
  end function H_dist_hot

  real (kind=8) function H_dist_hot_2(seed, nparam, hparam)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine draws randomly a number according to the hot belt H_r
! distribution, represented by an exponentially tapered exponential,
! with parameters I've fitted on the OSSOS hot belt data, in range 6-8.5.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : October 2021
! Version 2 : December 2021- updated parameter values.
! Version 3 : January 2022 - forcing slope at small sizes.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     nparam: Number of parameters (I4)
!     hparam: Parameters for asymptotic slope(s) (n*R8)
!             hparam(1): start of asymptote
!             hparam(2): contrast at start of asymptote
!             hparam(3): slope of asymptote
!             hparam(4): end of asymptote
!
! OUTPUT
!     H_dist_hot_2 : Random value of H (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) seed
    implicit none

    integer (kind=4), parameter :: np = 16384
    integer (kind=4), intent(inout) :: seed
    integer (kind=4), intent(in) :: nparam
    real (kind=8), intent(in) :: hparam(4)
    integer (kind=4) :: i
    real (kind=8) :: params(4), random, h_min, h_max
    real (kind=8), save :: proba(0:np), htab(0:np)
    real (kind=8) :: n1, n2, n3, sl1, sl2, sl3, c1, cb, c, h1, h2, h3
    logical, save :: first

    data first /.true./
    data params /-2.465d0, 7.114d0, 0.666d0, 0.785d0/
    data h_min /-1.d0/
    data sl1 /0.13d0/, sl2 /0.6d0/
    data h1 /3.2d0/, h2 /6.d0/, n1 /3.d0/

    if (first) then
       h_max = hparam(nparam)
       htab(0) = h_min
       proba(0) = 1.d-10
       h3 = hparam(1)
       n2 = Variably_tapered(h2, params)
       n3 = Variably_tapered(h3, params)
       c = hparam(2)
       sl3 = hparam(3)
! The normalisation is done with the exponentially tapered exponential,
! as fitted on the OSSOS hot component. This determines the normalisation
! of the exponential between H = h2 and H = h1
! N(<H) = n2*10**(sl2*(H-h2))
! Then, there is an excess divot at h1. See
! [[file:///home/petit/Research/OSSOS/tes/OSSOSpapers/Papers/GlobalLuminosityFunction/CumDiffDistributions.org]]
! for the appropriate formula.
       cb = n1*sl1*log(10.d0)
       c1 = (n2-n1)*sl2/(n1*sl1*(10.d0**(sl2*(h2-h1))-1.d0))
       do i = 1, np
          htab(i) = h_min + dfloat(i)*(h_max-h_min)/dfloat(np)
          if (htab(i) .lt. h1) then
             proba(i) = n1*10.d0**(sl1*(htab(i)-h1))
          else if (htab(i) .lt. h2) then
             proba(i) = n1 + c1*cb &
                  *(10.d0**(sl2*(htab(i)-h1))-1.d0)/(sl2*log(10.d0))
          else if (htab(i) .lt. h3) then
             proba(i) = Variably_tapered(htab(i), params)
          else
             proba(i) = n3 + n3*c*(10.d0**(sl3*(htab(i)-h3)) - 1.d0)
          end if
       end do
       do i = 0, np
          proba(i) = proba(i)/proba(np)
       end do
       first = .false.
    end if

    random = ran_3(seed)
    H_dist_hot_2 = interp(proba, htab, random, np+1)

    return
  end function H_dist_hot_2

  subroutine H_cum_hot_3(nparam, hparam, np, hs, dist)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine computes the cumulative distribution of the hot belt H_r,
! represented by an exponentially tapered exponential, with parameters
! fitted on the OSSOS cold belt data, then scaled to hot.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : June 2024 - From H_dist_hot_3
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     hparam: Parameters for asymptotic slope(s) (n*R8)
!             hparam(1): start of asymptote
!             hparam(2): contrast at start of asymptote
!             hparam(3): slope of asymptote
!             hparam(4): end of asymptote
!     np    : Number of points to return
!     hs    : Values of H at which we want the distribution
!
! OUTPUT
!     hs(0) : H_min
!     dist  : Values of the cumulative distribution
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: hparam
!f2py intent(in) np
!f2py intent(in, out) hs
!f2py intent(out) dist
    implicit none

    integer (kind=4), intent(in) :: nparam, np
    real (kind=8), intent(in) :: hparam(4)
    real (kind=8), intent(inout) :: hs(0:np)
    real (kind=8), intent(out) :: dist(0:np)
    integer (kind=4) :: i
    real (kind=8), save :: params(4), h_min, h_max
    real (kind=8), save :: n1, n2, n3, sl1, sl2, sl3, c1, cb, c, h1, h2, h3, scale

    data params /-2.6d0, 8.1d0, 0.666d0, 0.42d0/
    data h_min /-1.d0/
    data sl1 /0.13d0/, sl2 /0.6d0/
    data h1 /3.2d0/, h2 /6.d0/, n1 /3.d0/
    data scale /2.2d0/

    h_max = hparam(nparam)
    h3 = hparam(1)
    n2 = scale*Variably_tapered(h2, params)
    n3 = scale*Variably_tapered(h3, params)
    c = hparam(2)
    sl3 = hparam(3)
! The normalisation is done with the exponentially tapered exponential,
! as fitted on the OSSOS cold component, then scaled by 2. This
! determines the normalisation of the exponential between H = h2
! and H = h1
! N(<H) = n2*10**(sl2*(H-h2))
! Then, there is an excess divot at h1. See
! [[file:///home/petit/Research/OSSOS/tes/OSSOSpapers/Papers/GlobalLuminosityFunction/CumDiffDistributions.org]]
! for the appropriate formula.
    cb = n1*sl1*log(10.d0)
    c1 = (n2-n1)*sl2/(n1*sl1*(10.d0**(sl2*(h2-h1))-1.d0))
    hs(0) = h_min
    dist(0) = 1.d-10
    do i = 1, np
       hs(i) = h_min + dfloat(i)*(h_max-h_min)/dfloat(np)
       if (hs(i) .lt. h1) then
          dist(i) = n1*10.d0**(sl1*(hs(i)-h1))
       else if (hs(i) .lt. h2) then
          dist(i) = n1 + c1*cb &
               *(10.d0**(sl2*(hs(i)-h1))-1.d0)/(sl2*log(10.d0))
       else if (hs(i) .lt. h3) then
          dist(i) = scale*Variably_tapered(hs(i), params)
       else
          dist(i) = n3 + n3*c*(10.d0**(sl3*(hs(i)-h3)) - 1.d0)
       end if
    end do

    return
  end subroutine H_cum_hot_3

  real (kind=8) function H_draw_hot_3(seed, nparam, hparam)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine draws randomly a number according to the hot belt H_r
! distribution, represented by a cumulative exponentially tapered exponential,
! with parameters fitted on the OSSOS cold belt data, then scaled to
! hot.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : June 2024 - From H_dist_hot_3
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     nparam: Number of parameters (I4)
!     hparam: Parameters for asymptotic slope(s) (n*R8)
!             hparam(1): start of asymptote
!             hparam(2): contrast at start of asymptote
!             hparam(3): slope of asymptote
!             hparam(4): end of asymptote
!
! OUTPUT
!     H_draw_hot_3 : Random value of H (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) seed
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: hparam
    implicit none

    integer (kind=4), parameter :: np = 16384
    integer (kind=4), intent(inout) :: seed
    integer (kind=4), intent(in) :: nparam
    real (kind=8), intent(in) :: hparam(4)
    integer (kind=4) :: i
    real (kind=8), save :: random, h_min, h_max
    real (kind=8), save :: proba(0:np), htab(0:np)
    logical, save :: first

    data first /.true./
    data h_min /-1.d0/

    if (first) then
       h_max = hparam(nparam)
       htab(0) = h_min
       do i = 1, np
          htab(i) = h_min + dfloat(i)*(h_max-h_min)/dfloat(np)
       end do
       call H_cum_hot_3(nparam, hparam, np, htab, proba)
       do i = 0, np
          proba(i) = proba(i)/proba(np)
       end do
       first = .false.
    end if

    random = ran_3(seed)
    H_draw_hot_3 = interp(proba, htab, random, np+1)

    return
  end function H_draw_hot_3

  subroutine H_diff_hot_4(nparam, hparam, np, hs, dist)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine computes the differential distribution of the hot belt H_r,
! represented by an exponentially tapered exponential, with parameters
! fitted on the OSSOS cold belt data, then scaled to hot.
! This version provides a continuous differential function, except for the divot.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : June 2024 - From H_dist_hot_4
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     hparam: Parameters for asymptotic slope(s) (n*R8)
!             hparam(1): start of asymptote
!             hparam(2): contrast at start of asymptote
!             hparam(3): slope of asymptote
!             hparam(4): end of asymptote
!     np    : Number of points to return
!     hs    : Values of H at which we want the distribution
!
! OUTPUT
!     hs(0) : H_min
!     dist  : Values of the cumulative distribution
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: hparam
!f2py intent(in) np
!f2py intent(in, out) hs
!f2py intent(out) dist
    implicit none

    integer (kind=4), intent(in) :: nparam, np
    real (kind=8), intent(in) :: hparam(4)
    real (kind=8), intent(inout) :: hs(0:np)
    real (kind=8), intent(out) :: dist(0:np)
    integer (kind=4) :: i
    real (kind=8), save :: params(4), h_min, h_max
    real (kind=8), save :: a1, a2, a3, sl1, sl2, sl3, h1, h2, h3, scale

    data params /-2.6d0, 8.1d0, 0.666d0, 0.42d0/
    data h_min /-1.d0/
    data sl1 /0.13d0/, sl2 /0.6d0/
    data h1 /3.2d0/, h2 /6.d0/
    data scale /2.2d0/

    h_max = hparam(nparam)
    h3 = hparam(1)
    a2 = scale*Variably_tapered_diff(h2, params)
    a3 = hparam(2)*scale*Variably_tapered_diff(h3, params)
    sl3 = hparam(3)
    a1 = a2*10.0d0**(sl2*(h1-h2))
! The normalisation is done with the exponentially tapered exponential,
! as fitted on the OSSOS cold component, then scaled by 2. This
! determines the normalisation of the exponential between H = h2
! and H = h1
! N(<H) = n2*10**(sl2*(H-h2))
! Then, there is an excess divot at h1. See
! [[file:///home/petit/Research/OSSOS/tes/OSSOSpapers/Papers/GlobalLuminosityFunction/CumDiffDistributions.org]]
! for the appropriate formula.
    do i = 0, np
       if (hs(i) .le. h_min) then
          dist(i) = 0.0d0
       else if (hs(i) .le. h1) then
          dist(i) = a1*10.0d0**(sl1*(hs(i)-h1))
       else if (hs(i) .le. h2) then
          dist(i) = a2*10.0d0**(sl2*(hs(i)-h2))
       else if (hs(i) .le. h3) then
          dist(i) = scale*Variably_tapered_diff(hs(i), params)
       else if (hs(i) .le. h_max) then
          dist(i) = a3*10.0d0**(sl3*(hs(i)-h3))
       else
          dist(i) = 0.0d0
       end if
    end do
    
    return
  end subroutine H_diff_hot_4

  real (kind=8) function H_draw_hot_4(seed, nparam, hparam)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine draws randomly a number according to the hot belt H_r
! distribution, represented by a differntial exponentially tapered exponential,
! with parameters fitted on the OSSOS cold belt data, then scaled to
! hot.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : June 2024 - From H_dist_hot_4
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     nparam: Number of parameters (I4)
!     hparam: Parameters for asymptotic slope(s) (n*R8)
!             hparam(1): start of asymptote
!             hparam(2): contrast at start of asymptote
!             hparam(3): slope of asymptote
!             hparam(4): end of asymptote
!
! OUTPUT
!     H_draw_hot_4 : Random value of H (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) seed
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: hparam
    implicit none

    integer (kind=4), parameter :: np = 16384
    integer (kind=4), intent(inout) :: seed
    integer (kind=4), intent(in) :: nparam
    real (kind=8), intent(in) :: hparam(*)
    integer (kind=4) :: i
    real (kind=8), save :: random, h_min, h_max
    real (kind=8), save :: proba(0:np), htab(0:np)
    logical, save :: first

    data first /.true./
    data h_min /-1.d0/

    if (first) then
       h_max = hparam(nparam)
       htab(0) = h_min
       do i = 1, np
          htab(i) = h_min + dfloat(i)*(h_max-h_min)/dfloat(np)
       end do
       call H_diff_hot_4(nparam, hparam, np, htab, proba)
       proba(0) = 0.0d0
       do i = 1, np
          proba(i) = (proba(i)+proba(i-1))*(htab(i)-htab(i-1))/2.0d0 + proba(i-1)
       end do
       proba(0) = 1.0d-10
       do i = 0, np
          proba(i) = proba(i)/proba(np)
       end do
       first = .false.
    end if

    random = ran_3(seed)
    H_draw_hot_4 = interp(proba, htab, random, np+1)

  end function H_draw_hot_4

  subroutine H_diff_hot_5(nparam, hparam, np, hs, dist)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine computes the differential distribution of the hot belt H_r,
! represented by an exponentially tapered exponential, with parameters
! fitted on the OSSOS cold belt data, then scaled to hot.
! This version provides a continuous differential function, except for the divot.
!
! This version has modified parameters to obtain a cumulative distribution
! that looks like, and has the same normalisation as H_dist_hot_3.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : June 2024 - From H_dist_hot_4
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     hparam: Parameters for asymptotic slope(s) (n*R8)
!             hparam(1): start of asymptote
!             hparam(2): contrast at start of asymptote
!             hparam(3): slope of asymptote
!             hparam(4): end of asymptote
!     np    : Number of points to return
!     hs    : Values of H at which we want the distribution
!
! OUTPUT
!     hs(0) : H_min
!     dist  : Values of the cumulative distribution
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: hparam
!f2py intent(in) np
!f2py intent(in, out) hs
!f2py intent(out) dist
    implicit none

    integer (kind=4), intent(in) :: nparam, np
    real (kind=8), intent(in) :: hparam(*)
    real (kind=8), intent(inout) :: hs(0:np)
    real (kind=8), intent(out) :: dist(0:np)
    integer (kind=4) :: i
    real (kind=8), save :: params(4), h_min, h_max
    real (kind=8), save :: a1, a2, a3, sl1, sl2, sl3, h1, h2, h3, scale

    data params /-2.6d0, 8.1d0, 0.666d0, 0.42d0/
    data h_min /-1.d0/
    data sl1 /0.13d0/, sl2 /0.7d0/
    data h1 /2.5d0/, h2 /5.8d0/
    data scale /2.2d0/

    h_max = hparam(nparam)
    h3 = hparam(1)
    a2 = scale*Variably_tapered_diff(h2, params)
    a3 = hparam(2)*scale*Variably_tapered_diff(h3, params)
    sl3 = hparam(3)
    a1 = a2*10.0d0**(sl2*(h1-h2))
! The normalisation is done with the exponentially tapered exponential,
! as fitted on the OSSOS cold component, then scaled by 2. This
! determines the normalisation of the exponential between H = h2
! and H = h1
! N(<H) = n2*10**(sl2*(H-h2))
! Then, there is an excess divot at h1. See
! [[file:///home/petit/Research/OSSOS/tes/OSSOSpapers/Papers/GlobalLuminosityFunction/CumDiffDistributions.org]]
    ! for the appropriate formula.
    do i = 0, np
       if (hs(i) .le. h_min) then
          dist(i) = 0.0d0
       else if (hs(i) .le. h1) then
          dist(i) = a1*10.0d0**(sl1*(hs(i)-h1))
       else if (hs(i) .le. h2) then
          dist(i) = a2*10.0d0**(sl2*(hs(i)-h2))
       else if (hs(i) .le. h3) then
          dist(i) = scale*Variably_tapered_diff(hs(i), params)
       else if (hs(i) .le. h_max) then
          dist(i) = a3*10.0d0**(sl3*(hs(i)-h3))
       else
          dist(i) = 0.0d0
       end if
    end do
    
    return
  end subroutine H_diff_hot_5

  real (kind=8) function H_draw_hot_5(seed, nparam, hparam)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine draws randomly a number according to the hot belt H_r
! distribution, represented by a differntial exponentially tapered exponential,
! with parameters fitted on the OSSOS cold belt data, then scaled to
! hot.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : June 2024 - From H_dist_hot_4
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     nparam: Number of parameters (I4)
!     hparam: Parameters for asymptotic slope(s) (n*R8)
!             hparam(1): start of asymptote
!             hparam(2): contrast at start of asymptote
!             hparam(3): slope of asymptote
!             hparam(4): end of asymptote
!
! OUTPUT
!     H_draw_hot_5 : Random value of H (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) seed
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: hparam
    implicit none

    integer (kind=4), parameter :: np = 16384
    integer (kind=4), intent(inout) :: seed
    integer (kind=4), intent(in) :: nparam
    real (kind=8), intent(in) :: hparam(*)
    integer (kind=4) :: i
    real (kind=8), save :: random, h_min, h_max
    real (kind=8), save :: proba(0:np), htab(0:np)
    logical, save :: first

    data first /.true./
    data h_min /-1.d0/

    if (first) then
       h_max = hparam(nparam)
       htab(0) = h_min
       do i = 1, np
          htab(i) = h_min + dfloat(i)*(h_max-h_min)/dfloat(np)
       end do
       call H_diff_hot_5(nparam, hparam, np, htab, proba)
       proba(0) = 0.0d0
       do i = 1, np
          proba(i) = (proba(i)+proba(i-1))*(htab(i)-htab(i-1))/2.0d0 + proba(i-1)
       end do
       proba(0) = 1.0d-10
       do i = 0, np
          proba(i) = proba(i)/proba(np)
       end do
       first = .false.
    end if

    random = ran_3(seed)
    H_draw_hot_5 = interp(proba, htab, random, np+1)

  end function H_draw_hot_5

  subroutine H_diff_hot_6(nparam, hparam, np, hs, dist)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine computes the differential distribution of the hot belt H_r,
! represented by an exponentially tapered exponential, with parameters
! fitted on the OSSOS cold belt data, then scaled to hot.
! This version provides a continuous differential function, except for the divot.
!
! This version has modified parameters to obtain a cumulative distribution
! that looks like, and has the same normalisation as H_dist_hot_3.
!
! Here, there is a simple divot at given mag, and then the distribution
! continues according to the expoential taper, simply rescaled by the contrast.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : June 2024 - From H_dist_hot_5
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     hparam: Parameters for asymptotic slope(s) (n*R8)
!             hparam(1): start of asymptote
!             hparam(2): contrast at divot
!             hparam(3): end of distribution
!     np    : Number of points to return
!     hs    : Values of H at which we want the distribution
!
! OUTPUT
!     hs(0) : H_min
!     dist  : Values of the cumulative distribution
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: hparam
!f2py intent(in) np
!f2py intent(in, out) hs
!f2py intent(out) dist
    implicit none

    integer (kind=4), intent(in) :: nparam, np
    real (kind=8), intent(in) :: hparam(*)
    real (kind=8), intent(inout) :: hs(0:np)
    real (kind=8), intent(out) :: dist(0:np)
    integer (kind=4) :: i
    real (kind=8), save :: params(4), h_min, h_max
    real (kind=8), save :: a1, a2, a3, sl1, sl2, sl3, h1, h2, h3, scale

    data params /-2.6d0, 8.1d0, 0.666d0, 0.42d0/
    data h_min /-1.d0/
    data sl1 /0.13d0/, sl2 /0.7d0/
    data h1 /2.5d0/, h2 /5.8d0/
    data scale /2.2d0/

!    print *, nparam
!    print *, hparam(1:nparam)
    h_max = hparam(nparam)
    h3 = hparam(1)
    a2 = scale*Variably_tapered_diff(h2, params)
    a1 = a2*10.0d0**(sl2*(h1-h2))
!    print *, a1, a2, a3, h3, h_max
!    print *, scale*Variably_tapered_diff(h3, params), hparam(2)*scale*Variably_tapered_diff(h3, params)
! The normalisation is done with the exponentially tapered exponential,
! as fitted on the OSSOS cold component, then scaled by 2. This
! determines the normalisation of the exponential between H = h2
! and H = h1
! N(<H) = n2*10**(sl2*(H-h2))
! Then, there is an excess divot at h1. See
! [[file:///home/petit/Research/OSSOS/tes/OSSOSpapers/Papers/GlobalLuminosityFunction/CumDiffDistributions.org]]
    ! for the appropriate formula.
    do i = 0, np
       if (hs(i) .le. h_min) then
          dist(i) = 0.0d0
       else if (hs(i) .le. h1) then
          dist(i) = a1*10.0d0**(sl1*(hs(i)-h1))
       else if (hs(i) .le. h2) then
          dist(i) = a2*10.0d0**(sl2*(hs(i)-h2))
       else if (hs(i) .le. h3) then
          dist(i) = scale*Variably_tapered_diff(hs(i), params)
       else if (hs(i) .le. h_max) then
          dist(i) = hparam(2)*scale*Variably_tapered_diff(hs(i), params)
       else
          dist(i) = 0.0d0
       end if
    end do
    
    return
  end subroutine H_diff_hot_6

  real (kind=8) function H_draw_hot_6(seed, nparam, hparam)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine draws randomly a number according to the hot belt H_r
! distribution, represented by a differntial exponentially tapered exponential,
! with parameters fitted on the OSSOS cold belt data, then scaled to
! hot.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : June 2024 - From H_dist_hot_5
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     nparam: Number of parameters (I4)
!     hparam: Parameters for asymptotic slope(s) (n*R8)
!             hparam(1): start of asymptote
!             hparam(2): contrast at divot
!             hparam(3): end of distribution
!
! OUTPUT
!     H_draw_hot_6 : Random value of H (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) seed
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: hparam
    implicit none

    integer (kind=4), parameter :: np = 16384
    integer (kind=4), intent(inout) :: seed
    integer (kind=4), intent(in) :: nparam
    real (kind=8), intent(in) :: hparam(*)
    integer (kind=4) :: i
    real (kind=8), save :: random, h_min, h_max
    real (kind=8), save :: proba(0:np), htab(0:np)
    logical, save :: first

    data first /.true./
    data h_min /-1.d0/

    if (first) then
       h_max = hparam(nparam)
       htab(0) = h_min
       do i = 1, np
          htab(i) = h_min + dfloat(i)*(h_max-h_min)/dfloat(np)
       end do
       call H_diff_hot_6(nparam, hparam, np, htab, proba)
       proba(0) = 0.0d0
       do i = 1, np
          proba(i) = (proba(i)+proba(i-1))*(htab(i)-htab(i-1))/2.0d0 + proba(i-1)
       end do
       proba(0) = 1.0d-10
       do i = 0, np
          proba(i) = proba(i)/proba(np)
       end do
       first = .false.
    end if

    random = ran_3(seed)
    H_draw_hot_6 = interp(proba, htab, random, np+1)

  end function H_draw_hot_6

  real (kind=8) function offgau (nparam, param, inc)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the unnormalized inclination "probability"
! density as a non-zero centered gaussian.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : August 2014
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     param : Parameters (n*R8)
!     inc   : Inclination [rad] (R8)
!
! OUPUT
!     offgau: Value of the probability (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: param
!f2py intent(in) inc
    implicit none

    integer (kind=4), intent(in) :: nparam
    real (kind=8), intent(in) :: param(*), inc
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi
    real (kind=8) :: fe, s1, angle, t0, t1, t3

    if (nparam .ne. 2) stop
    t0 = param(1)
    s1 = param(2)
    t1 = 2.*s1**2
    angle = mod(inc, TwoPi)
    if (angle .gt. Pi) angle = angle - TwoPi
    t3 = -(t0-angle)**2
    if (t3 .lt. -300.d0*t1) then
       fe = 0.d0
    else
       fe = exp(t3/t1)
    end if
    offgau = dsin(angle)*fe

    return
  end function offgau

  real (kind=8) function cold_low_a_inc (nparam, param, inc)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the unnormalized inclination "probability"
! density for cold objects with a < 44.4
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : December 2019
! Version 2 : December 2019, 20th
! Version 3 : December 2019, 28th
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     param : Parameters (n*R8)
!     inc   : Inclination [rad] (R8)
!
! OUPUT
!     cold_low_a_inc: Value of the probability (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: param
!f2py intent(in) inc
    implicit none

    integer (kind=4), intent(in) :: nparam
    real (kind=8), intent(in) :: param(*), inc
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, &
         TwoPi = 2.0d0*Pi, drad = Pi/180.0d0
    real (kind=8) :: angle

    angle = mod(inc, TwoPi)
    if (angle .gt. Pi) angle = angle - TwoPi
    if (angle .lt. 0.8d0*drad) then
       cold_low_a_inc = angle/(0.8*drad)
    else if (angle .lt. 2.6d0*drad) then
       cold_low_a_inc = 1.d0
    else if (angle .le. 4.3d0*drad) then
       cold_low_a_inc = 0.7d0*(4.3d0*drad-angle)/(1.7d0*drad) + 0.3d0
    else if (angle .le. 4.7d0*drad) then
       cold_low_a_inc = 0.3d0*(4.7d0*drad - angle)/(0.4d0*drad)
    else
       cold_low_a_inc = 0.d0
    end if
      
    return
  end function cold_low_a_inc

  real (kind=8) function cold_low_a_inc_2 (nparam, param, inc)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the unnormalized inclination "probability"
! density for cold objects with a < 44.4
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : December 2019
! Version 2 : December 2019, 20th
! Version 3 : December 2019, 28th
! Version 4 : December 2021, 6th
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     param : Parameters (n*R8)
!     inc   : Inclination [rad] (R8)
!
! OUPUT
!     cold_low_a_inc_2: Value of the probability (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: param
!f2py intent(in) inc
    implicit none

    integer (kind=4), intent(in) :: nparam
    real (kind=8), intent(in) :: param(*), inc
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, &
         TwoPi = 2.0d0*Pi, drad = Pi/180.0d0
    real (kind=8) :: angle
    real (kind=8), parameter :: a1 = 0.3d0, a2 = 4.05d0, a3 = 4.5d0

    angle = mod(inc, TwoPi)
    if (angle .gt. Pi) angle = angle - TwoPi
    if (angle .lt. a1*drad) then
       cold_low_a_inc_2 = angle/(a1*drad)
    else if (angle .lt. a2*drad) then
       cold_low_a_inc_2 = 1.d0
    else if (angle .le. a3*drad) then
       cold_low_a_inc_2 = (a3*drad - angle)/((a3-a2)*drad)
    else
       cold_low_a_inc_2 = 0.d0
    end if
      
    return
  end function cold_low_a_inc_2

  real (kind=8) function kernel_inc (nparam, param, inc)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the unnormalized inclination "probability"
! density for cold objects in the kernel region.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : December 2021, 7th
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     param : Parameters (n*R8)
!     inc   : Inclination [rad] (R8)
!
! OUPUT
!     kernel_inc: Value of the probability (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: param
!f2py intent(in) inc
    implicit none

    integer (kind=4), intent(in) :: nparam
    real (kind=8), intent(in) :: param(*), inc
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, &
         TwoPi = 2.0d0*Pi, drad = Pi/180.0d0
    real (kind=8), save :: angle, sum_la, da, lli, lhi
    real (kind=8) :: fkrla, fkla, fli
    integer (kind=4), save :: i, nstep
    logical, save :: first
    real (kind=8), parameter :: a1 = 0.9d0, a2 = 2.0d0, delta = 0.8d0

    common /com_kernel/ fkrla, fkla, fli

    data first /.true./, nstep /5000/

    if (first) then
       sum_la = 0.d0
       da = 5.d0*drad/dble(nstep)
       do i = 1, nstep
          angle = i*da
          sum_la = sum_la + cold_low_a_inc_2(nparam, param, angle)
       end do
       sum_la = sum_la*da*fkla/(fkrla-fkla)
       lli = fli*sum_la/(delta*drad)
       lhi = (1.d0-fli)*sum_la/(delta*drad)
    end if

    angle = mod(inc, TwoPi)
    if (angle .gt. Pi) angle = angle - TwoPi
    kernel_inc = cold_low_a_inc_2(nparam, param, angle)
    if ((angle .ge. a1*drad) .and. (angle .le. (a1+delta)*drad)) then
       kernel_inc = kernel_inc + lli
    else if ((angle .ge. a2*drad).and.(angle .le. (a2+delta)*drad)) then
       kernel_inc = kernel_inc + lhi
    end if
      
    return
  end function kernel_inc

  real (kind=8) function cold_inc_sq (nparam, param, inc)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the unnormalized inclination "probability"
! density for cold objects.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : April 2023, 7th
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     param : Parameters (n*R8)
!     inc   : Inclination [rad] (R8)
!
! OUPUT
!     cold_inc_sq: Value of the probability (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: param
!f2py intent(in) inc
    implicit none

    integer (kind=4), intent(in) :: nparam
    real (kind=8), intent(in) :: param(*), inc
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, &
         TwoPi = 2.0d0*Pi, drad = Pi/180.0d0
    real (kind=8), save :: angle, a0, p
    logical, save :: first

    data a0 /5.0d0/, p /1.0d0/
    data first /.true./

    if (first) then
       a0 = a0*drad
       if (nparam .ge. 1) then
          a0 = param(1)
       end if
       if (nparam .ge. 2) then
          p = param(2)
       end if
       first = .false.
    end if

    angle = mod(inc, TwoPi)
    if ((angle .gt. 0.0d0) .and. (angle .lt. a0)) then
       cold_inc_sq = angle**p*(a0**p - angle**p)
    else
       cold_inc_sq = 0.0d0
    end if

    return
  end function cold_inc_sq

  real (kind=8) function warm_inc_sq (nparam, param, inc)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the unnormalized inclination "probability"
! density for warm objects.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : April 2023, 7th
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     param : Parameters (n*R8)
!     inc   : Inclination [rad] (R8)
!
! OUPUT
!     warm_inc_sq: Value of the probability (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: param
!f2py intent(in) inc
    implicit none

    integer (kind=4), intent(in) :: nparam
    real (kind=8), intent(in) :: param(*), inc
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, &
         TwoPi = 2.0d0*Pi, drad = Pi/180.0d0
    real (kind=8), save :: angle, a0, p
    logical, save :: first

    data a0 /12.0d0/, p /1.0d0/
    data first /.true./

    if (first) then
       a0 = a0*drad
       if (nparam .ge. 1) then
          a0 = param(1)
       end if
       if (nparam .ge. 2) then
          p = param(2)
       end if
       first = .false.
    end if

    angle = mod(inc, TwoPi)
    if ((angle .gt. 0.0d0) .and. (angle .lt. a0)) then
       warm_inc_sq = angle**p*(a0**p - angle**p)
    else
       warm_inc_sq = 0.0d0
    end if

    return
  end function warm_inc_sq

  real (kind=8) function warm_inc_sk_sq (nparam, param, inc)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the unnormalized inclination "probability"
! density for warm objects.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : April 2023, 7th
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     param : Parameters (n*R8)
!     inc   : Inclination [rad] (R8)
!
! OUPUT
!     warm_inc_sk_sq: Value of the probability (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: param
!f2py intent(in) inc
    implicit none

    integer (kind=4), intent(in) :: nparam
    real (kind=8), intent(in) :: param(*), inc
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, &
         TwoPi = 2.0d0*Pi, drad = Pi/180.0d0
    real (kind=8), save :: angle, a0, p
    logical, save :: first

    data a0 /10.0d0/, p /1.0d0/
    data first /.true./

    if (first) then
       a0 = a0*drad
       if (nparam .ge. 1) then
          a0 = param(1)
       end if
       if (nparam .ge. 2) then
          p = param(2)
       end if
       first = .false.
    end if

    angle = mod(inc, TwoPi)
    if ((angle .gt. 0.0d0) .and. (angle .lt. a0)) then
       warm_inc_sk_sq = sin(angle*Pi/(a0*2.0d0))**p*(1.0d0 - sin(angle*Pi/(a0*2.0d0))**p)
    else
       warm_inc_sk_sq = 0.0d0
    end if

    return
  end function warm_inc_sk_sq

  real (kind=8) function cold_high_a_inc (nparam, param, inc)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the unnormalized inclination "probability"
! density for cold objects with a < 44.4
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : December 2019
! Version 2 : December 2019, 20th
! Version 3 : December 2019, 28th
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     param : Parameters (n*R8)
!     inc   : Inclination [rad] (R8)
!
! OUPUT
!     cold_high_a_inc: Value of the probability (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: param
!f2py intent(in) inc
    implicit none

    integer (kind=4), intent(in) :: nparam
    real (kind=8), intent(in) :: param(*), inc
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, &
         TwoPi = 2.0d0*Pi, drad = Pi/180.0d0
    real (kind=8) :: angle
    real (kind=8), parameter :: a1 = 1.2d0, a2 = 4.0d0, a3 = 6.0d0, a4 = 12.0d0

    angle = mod(inc, TwoPi)
    if (angle .gt. Pi) angle = angle - TwoPi
    if (angle .lt. a1*drad) then
       cold_high_a_inc = angle/(a1*drad)
    else if (angle .le. a2*drad) then
       cold_high_a_inc = 1.d0
    else if (angle .le. a3*drad) then
       cold_high_a_inc = 0.3d0*(a3*drad-angle)/((a3-a2)*drad) + 0.7d0
    else if (angle .le. a4*drad) then
       cold_high_a_inc = 0.7d0*(a4*drad - angle)/((a4-a3)*drad)
    else
       cold_high_a_inc = 0.d0
    end if

    return
  end function cold_high_a_inc

  real (kind=8) function hot_inc (nparam, param, inc)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the unnormalized inclination "probability"
! density for cold objects with a < 44.4
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : December 2019
! Version 2 : December 2019, 20th
! Version 3 : December 2019, 28th
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     param : Parameters (n*R8)
!     inc   : Inclination [rad] (R8)
!
! OUPUT
!     hot_inc: Value of the probability (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: param
!f2py intent(in) inc
    implicit none

    integer (kind=4), intent(in) :: nparam
    real (kind=8), intent(in) :: param(*), inc
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, &
         TwoPi = 2.0d0*Pi, drad = Pi/180.0d0, PiThird = Pi/3.0d0
    real (kind=8) :: angle
    real (kind=8), save :: a0, a1, a2, a3
    logical, save :: first

    data a0 /5.0d0/, a1 /7.0d0/, a2 /26.0d0/, a3 /46.0d0/
    data first /.true./

    if (first) then
       a0 = a0*drad
       a1 = a1*drad
       a2 = a2*drad
       a3 = a3*drad
       if (nparam .ge. 1) then
          a0 = param(1)
       end if
       if (nparam .ge. 2) then
          a1 = param(2)
       end if
       if (nparam .ge. 3) then
          a2 = param(3)
       end if
       if (nparam .ge. 4) then
          a3 = param(4)
       end if
       first = .false.
    end if

    angle = mod(inc, TwoPi)
    if (angle .gt. Pi) angle = angle - TwoPi
    if (angle .lt. a0) then
       hot_inc = 0.d0
    else if (angle .lt. a1) then
       hot_inc = 0.5d0*(angle-a0)/(a1-a0)
    else if (angle .le. a2) then
       hot_inc = 0.5d0
    else if (angle .le. a3) then
       hot_inc = 1.d0 - cos((a3-angle)/(a3-a2)*PiThird)
    else
       hot_inc = 0.d0
    end if
 
    return
  end function hot_inc

  real (kind=8) function onecomp (nparam, param, inc)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the unnormalized inclination "probability"
! density of Brown.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : February 2007
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     param : Parameters (n*R8)
!     inc   : Inclination [rad] (R8)
!
! OUPUT
!     onecomp: Value of the probability (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: param
!f2py intent(in) inc
    implicit none

    integer (kind=4), intent(in) :: nparam
    real (kind=8), intent(in) :: param(*), inc
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi
    real (kind=8) :: fe, s1, angle, t1, t3

    if (nparam .ne. 1) stop
    s1 = param(1)
    t1 = 2.*s1**2
    angle = mod(inc, TwoPi)
    if (angle .gt. Pi) angle = angle - TwoPi
    t3 = -angle**2
    if (t3 .lt. -300.d0*t1) then
       fe = 0.d0
    else
       fe = exp(t3/t1)
    end if
    onecomp = dsin(angle)*fe

    return
  end function onecomp

  real (kind=8) function onecompjmp (nparam, param, inc)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the unnormalized inclination "probability" density of
! Brown distribution function, modified to include sin(i)**2 instead of sin(i).
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : April 2020
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     param : Parameters (n*R8)
!     inc   : Inclination [rad] (R8)
!
! OUPUT
!     onecompjmp: Value of the probability (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: param
!f2py intent(in) inc
    implicit none

    integer (kind=4), intent(in) :: nparam
    real (kind=8), intent(in) :: param(*), inc
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, TwoPi = 2.0d0*Pi
    real (kind=8) :: fe, s1, angle, t1, t3

    if (nparam .ne. 1) stop
    s1 = param(1)
    t1 = 2.*s1**2
    angle = mod(inc, TwoPi)
    if (angle .gt. Pi) angle = angle - TwoPi
    t3 = -angle**2
    if (t3 .lt. -300.d0*t1) then
       fe = 0.d0
    else
       fe = exp(t3/t1)
    end if
    onecompjmp = dsin(angle)**2*fe

    return
  end function onecompjmp

  real (kind=8) function interp (x, y, val, n)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This function linearly interpolates the function y(x) at value x=val.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : October 2006
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     x     : Abscissa of the function, sorted in ascending order (n*R8)
!     y     : Values of the function (n*R8)
!     val   : Value of x at which to interpolate (R8)
!     n     : Size of x and y arrays (I4)
!
! OUTPUT
!     interp: Interpolated value (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) n
!f2py intent(in), depend(n) :: x
!f2py intent(in), depend(n) :: y
!f2py intent(in) val
    implicit none

    integer (kind=4), intent(in) ::  n
    real (kind=8), intent(in) :: x(*), y(*), val
    integer ::  ilo, ihi, i

    if (val .le. x(1)) then
       interp = y(1)
    else if (val .ge. x(n)) then
       interp = y(n)
    else
       ilo = 1
       ihi = n
1000   continue
       if (ihi - ilo .gt. 1) then
          i = (ihi + ilo)/2
          if (x(i) .lt. val) then
             ilo = i
          else if (x(i) .gt. val) then
             ihi = i
          else
             interp = y(i)
             return
          end if
          goto 1000
       end if
       interp = y(ilo) + (y(ihi) - y(ilo))*(val - x(ilo))/(x(ihi) - x(ilo))
    end if

    return
  end function interp

  real (kind=8) function size_dist_one (seed, h_params)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This function draws a number according to an exponential differential
! distribution specified by "h_params":
!
!   P(h) d_h = A \exp{(h*h_params(3))} d_h
!
! with h_params(1) <= h <= h_params(2)
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : April 2004
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : seed for the random number generator (I4)
!     h_params: parameters for the distribution (3*R8)
!
! OUTPUT
!     size_dist_one: Random value of H (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! Set of F2PY directives to create a Python module
!
!f2py intent(in,out) seed
!f2py intent(in) h_param
    implicit none

    integer (kind=4), intent(inout) :: seed
    real (kind=8), intent(in) :: h_params(*)
    real (kind=8) :: h0s10, h1s10, random, slope
!
! H-mag distribution 
!
! Functions have been triple-checked as of 2018-05-04. OK.
    slope = h_params(3)
    h0s10 = 10.d0**(slope*h_params(1))
    h1s10 = 10.d0**(slope*h_params(2))
    random=ran_3(seed)
    size_dist_one = log10( random*(h1s10 - h0s10) + h0s10 ) / slope

    return
  end function size_dist_one

  real (kind=8) function size_dist_two (seed, h_params)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This function draws a number according to a 2slope exponential
! differential distribution specified by "h_params":
!
!   P(h) d_h = A \exp{(h*h_params(4))} d_h
!
! with h_params(1) <= h <= h_params(2)
!
!   P(h) d_h = B \exp{(h*h_params(5))} d_h
!
! with h_params(2) <= h <= h_params(3)
!
! continuous at H1 = h_params(2)
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : April 2014
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : seed for the random number generator (I4)
!     h_params: parameters for the distribution (3*R8)
!
! OUTPUT
!     size_dist_two: Random value of H (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! Set of F2PY directives to create a Python module
!
!f2py intent(in,out) seed
!f2py intent(in) h_param
    implicit none

    integer (kind=4), intent(inout) :: seed
    real (kind=8), intent(in) :: h_params(5)
    real (kind=8) :: h0s1, h1s1, h1s2, h2s2, h1s12, random, sl1, sl2, xi1
!
! H-mag distribution
!
! Functions have been triple-checked as of 2018-05-04. OK.
    sl1 = h_params(4)
    sl2 = h_params(5)
    h0s1 = 10.d0**(sl1*h_params(1))
    h1s1 = 10.d0**(sl1*h_params(2))
    h1s2 = 10.d0**(sl2*h_params(2))
    h1s12 = h1s1/h1s2
    h2s2 = 10.d0**(sl2*h_params(3))
    xi1 = sl2*(h1s1-h0s1)/(sl1*h1s12*(h2s2-h1s2)+sl2*(h1s1-h0s1))
    random=ran_3(seed)
    if (random .le. xi1) then
       size_dist_two = log10(random/xi1*(h1s1 - h0s1) + h0s1)/sl1
    else
       size_dist_two = log10((random-xi1)/(1.d0-xi1)*(h2s2 - h1s2) &
            + h1s2)/sl2
    end if

    return
  end function size_dist_two

  real (kind=8) function size_dist_two_straight (seed, h_params)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This function draws a number according to a 2-slope exponential
! cumulative distribution specified by "h_params" (2 straight lines in
! semilog for cumulative):
!
!   P(<=h) = A 10^{(h*h_params(4))}
!
! with h_params(1) <= h <= h_params(2)
!
!   P(<=h) = B 10^{(h*h_params(5))}
!
! with h_params(2) < h <= h_params(3)
!
! continuous at H1 = h_params(2)
!   (A = B 10^{(h_params(5)-h_params(4))*h_params(2)})
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : May 2018
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : seed for the random number generator (I4)
!     h_params: parameters for the distribution (3*R8)
!
! OUTPUT
!     size_dist_two_straight: Random value of H (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! Set of F2PY directives to create a Python module
!
!f2py intent(in,out) seed
!f2py intent(in) h_param
    implicit none

    integer (kind=4), intent(inout) :: seed
    real (kind=8), intent(in) :: h_params(5)
    real (kind=8) :: h0s1, h1s1, h1s2, h2s2, h12s2, random, sl1, sl2
!
! H-mag distribution
!
! Functions have been triple-checked as of 2018-05-04. OK.
!
! In order to not have the ramp-up found in {\it size_dist_two}, I use
! an exponential with no lower limit on H. What sets the lower limit of
! the faint end is the selection of the correct slope at H_break (see
! test below). The test is based on the number of objects in each
! portion of the distribution.
!
! Since we want N(<H) = A 10^{\alpha_1 H} for H < H_b
! and           N(<H) = B 10^{\alpha_2 H} for H >= H_b
! and N(<H) continuous at H_b, we have N(<H_b) = A 10^{\alpha_1 H_b} =
! B 10^{\alpha_2 H_b}.
!
! The total number of objects is N(<H_max) = B 10^{\alpha_2 H_max}, and
! the number of objects bigger than H_b is N(<H_b) = B 10^{\alpha_2 H_b}. 
! So the fraction of objects bigger than H_b is 10^{\alpha_2 H_b} /
! 10^{\alpha_2 H_max} = h1s2 / h2s2 = h12s2.
!
! Therefore, if 'random' < h12s2, then we have a big object, and
! otherwise a small on. 'random' needs to be rescaled to go to 1.
    sl1 = h_params(4)
    sl2 = h_params(5)
    h0s1 = 10.d0**(sl1*h_params(1))
    h1s1 = 10.d0**(sl1*h_params(2))
    h1s2 = 10.d0**(sl2*h_params(2))
    h2s2 = 10.d0**(sl2*h_params(3))
    h12s2 = h1s2/h2s2
    random=ran_3(seed)
    if (random .le. h12s2) then
       size_dist_two_straight = log10(random*h1s1/h12s2)/sl1
    else
       size_dist_two_straight = log10(random*h2s2)/sl2
    end if

    return
  end function size_dist_two_straight

  real (kind=8) function size_dist_n_straight (seed, h_params, nn)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This function draws a number according to an n-slope exponential
! cumulative distribution specified by "h_params" (n straight lines in
! semilog for cumulative):
!
!   P(<=h) = A_k 10^{(h*\alpha_k)}
!
! with H_{k-1} < h <= H_k,
!
! with
!
! \alpha_k = h_params(k+n)
! H_k = h_params(k)
!
! for k in [1; n ]. h_params(0) is not present and implicit at -\infty.
! This corresponds to dropping the lower limit as we want straight lines.
!
! The function is continuous at H_k = h_params(k) for all k's:
!   A_k 10^{H_k*\alpha_k} = A_{k+1} 10^{H_k*\alpha_{k+1}}
!
! for k in [1; n-1].
!
! To avoid using allocatable arrays and dynamical allocation, I restrict
! the number of slopes to be <= 10.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : September 2019
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : seed for the random number generator (I4)
!     h_params: parameters for the distribution (3*R8)
!     nn    : number of different slopes (I4)
!
! OUTPUT
!     size_dist_n_straight: Random value of H (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! Set of F2PY directives to create a Python module
!
!f2py intent(in,out) seed
!f2py intent(in) h_param
!f2py intent(in) nn
    implicit none

    integer (kind=4), intent(in) :: nn
    integer (kind=4), intent(inout) :: seed
    real (kind=8), intent(in) :: h_params(*)
    integer (kind=4), save :: n, k
    real (kind=8) :: random, x(10), f(10)
    logical, save :: first

    data first /.true./

    if (nn .gt. 10) then
       if (first) then
          first = .false.
          print *, &
               'WARNING: number of breaks greater than allowed maximum'// &
               ' 10, using only the 10 first values.'
       end if
       n = 10
    else
       n = nn
    end if
!
! H-mag distribution
!
! In order to not have the ramp-up found in {\it size_dist_two}, I use
! an exponential with no lower limit on H. What sets the lower limit of
! the faint end is the selection of the correct slope at H_k (see
! test below). The test is based on the number of objects in each
! portion of the distribution.
!
! The total number of objects is N(<H_max) = N(<H_n) = A_n 10^{H_n \alpha_n},
! and the number of objects bigger than the next break H_{n-1} is N(
! <H_{n-1}) = A_n 10^{H_{n-1} \alpha_n}. More generaly, the number of
! objects bigger than H_k is
!
! N_k = A_k 10^{H_k \alpha_k} = A_{k+1} 10^{H_k \alpha_{k+1}}
!
! Let's define the fraction of objects bigger than H_k in the objects
! bigger than H_{k+1} as:
!
! x_k = N_{k-1} / N_k = 10^{(H_{k-1}-H_k) \alpha_k}
!
! for k in [2; n ], and x_1 = 0.
!
! The fraction of object bigger than B_k compared to the total
! number of objects is
!
! f_k = N_k / N_n = (N_k / N_{k+1}) (N_{k+1} / N_{k+2}) ... (N_{n-1} / N_n)
!                 = x_{k+1} x_{k+2} ... x_n
!
! We define x_k as described above, then we set
!
! f_n = 1.
! f_{k-1} = x_k f_k, for k in [2; n ]
!
! Therefore, if (f_{k-1} < random <= f_k), then we have an object in
! range ]H_{k-1}; H_k]. 'random' needs to be rescaled to go to 1.
    x(1) = 0.d0
    do k = 2, n
       x(k) = 10.d0**((h_params(k-1)-h_params(k))*h_params(k+n))
    end do
    f(n) = 1.d0
    do k = n, 2, -1
       f(k-1) = f(k)*x(k)
    end do
    random=ran_3(seed)
    do k = 1, n
       if (random .le. f(k)) then
          size_dist_n_straight = &
               log10(random/f(k)*10.d0**(h_params(k)*h_params(k+n))) &
               /h_params(k+n)
          return
       end if
    end do

  end function size_dist_n_straight

  real (kind=8) function scat_a (seed, np, p)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine draws randomly a number according to the scattering a
! distribution, represented by 2 logarythms with parameterss I can adjust.
! The first values have been eye-balled on the scattering a distribution.
!
! The cumulative probability of $a$ is
! \begin{eqnarray}
! P(a) & = & A_1 + \beta_1 Log{(a)} \qquad {\rm for} \qquad a < a_b \\
! P(a) & = & A_2 + \beta_2 Log{(a)} \qquad {\rm for} \qquad a \ge a_b
! \end{eqnarray}
! where $a_b$ and $P_b = P(a_b)$ define the break and are given/adjusted $a_b = 100$
! and $P_b = 0.8$. Next, we can either define the maximum value of $a$ where $P =
! 1$ and the minimum $a$ where $P$ is the minimumvalue corresponding to 1
! occurance. Similarly, one can give the values of $\beta_1$ and $\beta_2$.
! It makes more sense to give $a_{mini}$ and $a_{maxi}$ and deduce $\beta_1$ and $\beta_2$.
! \begin{eqnarray}
! \beta_2 & = & \frac{P_b - 1}{\Log{(a_b/a_{maxi})} \\
! \beta_1 & = & \frac{P_b}{\Log{(a_b/a_{mini})}
! \end{eqnarray}
! This will set the values of $A_1$ and $A_2$ with $a_b$ and $P_b$.
! \begin{eqnarray}
! A_1 & = & P_b - \beta_1 \Log{(a_b)} \\
! A_2 & = & P_b - \beta_2 \Log{(a_b)}
! \end{eqnarray}
! One first draws $P$ uniformly in [0; 1], then computes $a$ from
! \begin{eqnarray}
! a & = & 10^{\frac{P - A_1}{\beta_1}} \qquad {\rm for} \qquad P < P_b \\
! a & = & 10^{\frac{P-A_2}{\beta_2}} \qquad {\rm for} \qquad P \ge P_b
! \end{eqnarray}
!
! Actually, once debiased, it turns out that it could be a single slope
! distribution. I'll keep the 2-slope shape, but set the parameters so it's
! like a single slope distribution.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : August 2023
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     np    : Number of parameters (I4)
!     p     : Parameters for distribution (n*R8)
!             p(1): break point $a_p$
!             p(2): cumulative proba at break, $P_b$
!             p(3): minimum value of a
!             p(4): maximum value of a
!
! OUTPUT
!     scat_a : Random value of a (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) seed
!f2py intent(in) np
!f2py intent(in) p
    implicit none

    integer (kind=4), intent(inout) :: seed
    integer (kind=4), intent(in) :: np
    real (kind=8), intent(in) :: p(*)
    integer (kind=4) :: i
    real (kind=8), save :: ab, pb, amini, amaxi, b1, b2, A1, A2, lab
    real (kind=8) :: random
    logical, save :: first

    data ab /700.0d0/, pb /1.0d0/, amini /30.0d0/, amaxi /705.0d0/
    data first /.true./

    if (first) then
       if (np .ge. 1) then
          ab = p(1)
       end if
       if (np .ge. 2) then
          pb = p(2)
       end if
       if (np .ge. 3) then
          amini = p(3)
       end if
       if (np .ge. 4) then
          amaxi = p(4)
       end if
       lab = log10(ab)
       b1 = pb/(lab-log10(amini))
       A1 = pb - b1*lab
       b2 = (pb-1.0d0)/(lab-log10(amaxi))
       A2 = pb - b2*lab
!       print *, 'Scat_a'
!       print *, ab, pb, amini, amaxi
!       print *, b1, A1, b2, A2
!       print *, 10.0**((0.0d0-A1)/b1), 10.0**((0.5d0-A1)/b1), 10.0**((0.9999d0-A1)/b1)
       first = .false.
    end if

    random=ran_3(seed)
    if (random .le. pb) then
       scat_a = 10.0**((random-A1)/b1)
    else
       scat_a = 10.0**((random-A2)/b2)
    end if

    return
  end function scat_a

  real (kind=8) function scat_a_2 (seed, np, p)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine draws randomly a number according to the scattering a
! distribution, represented by a power-law of the shifted semimajor axis.
! The first parameter values have been eye-balled on the scattering a distribution.
!
! The cumulative probability of $a$ is
! \begin{displaymath}
! P(a) = \left(\frac{a - a_{min}}{a_{max} - a_{min}}\right)^\alpha
! \end{displaymath}
! where $a_{min}$ and $a_{max}$ are the smallest and largest semimajor axis
! and $\alpha$ is the index that will set the global shape.
! 
! Once a probablity $P(a)$ is drawn in [0; 1], one can easily obtain $a$ with:
! \begin{displaymath}
! a = a_{min}} + (a_{max} - a_{min})P(a)^{1/\alpha}
! \end{displaymath}
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : August 2023
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     np    : Number of parameters (I4)
!     p     : Parameters for distribution (n*R8)
!             p(1): minimum value of a
!             p(2): maximum value of a
!             p(3): exponential slope of distribution
!
! OUTPUT
!     scat_a_2: Random value of a (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) seed
!f2py intent(in) np
!f2py intent(in) p
    implicit none

    integer (kind=4), intent(inout) :: seed
    integer (kind=4), intent(in) :: np
    real (kind=8), intent(in) :: p(*)
    integer (kind=4) :: i
    real (kind=8), save :: amini, amaxi, alpha, al1
    real (kind=8) :: random
    logical, save :: first

    data amini /30.0d0/, amaxi /700.0d0/, alpha /0.5d0/
    data first /.true./

    if (first) then
       if (np .ge. 1) then
          amini = p(1)
       end if
       if (np .ge. 2) then
          amaxi = p(2)
       end if
       if (np .ge. 3) then
          alpha = p(3)
       end if
       al1 = 1.0d0/alpha
!       print *, 'Scat_a_2'
!       print *, amini, amaxi, alpha
!       print *, al1
       first = .false.
    end if

    random=ran_3(seed)
    scat_a_2 = amini + (amaxi - amini)*random**al1

    return
  end function scat_a_2

  real (kind=8) function scat_q (seed, np, p, a)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine draws randomly a number according to the scattering q
! distribution, described below.
! The first values have been eye-balled on the scattering q distribution.
!
! The extend of the $q$ distribution depends on $a$. So we first
! scale according to $a$. The minimum value of $q$ is $q_{min} = 10$. The maximum
! value of $q$ is $q_{max} \propto a^\alpha$. From the detectiosn, I determined $a_1
! = 43$ corresponds to $q_1 = 33$ and $a_2 = 425$ corresponds to $q_2 = 46$. This
! sets $\alpha$ (which could be adjusted) and then the proportionnality factor. now
! I define a variable that runs in [0; 1]:
! \begin{displaymath}
! \frac{xq - xq_0}{xq_1 - xq_0}= \frac{q-q_0}{q_1 \left(\frac{a}{a_1}\right)^\alpha - q_0}
! \end{displaymath}
! where $xq_0 = 0$, $xq_1 = 1$, $q_0 = q_{min}$. Now, the distribution of $xq$ looks
! like an exponential. So I can set $\Log{(P)} = a \times xq + b$ where $P$ is tyhe
! cumulative probability. We have $P(xq=1) = 1$, thus $b = -a$. In this way,
! $P$ varies from $10^a = A^{}$ to 1. So I introduce anothere variable $\Chi =
! \frac{P - A}{1 - A}$, or inversly $P = (1 - A) \Chi + A$.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : August 2023
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     np    : Number of parameters (I4)
!     p     : Parameters for distribution (n*R8)
!             p(1): exponential slope of $xq$ distribution
!             p(2): low a value of $a, q$ curve
!             p(3): low q value of $a, q$ curve
!             p(4): high a value of $a, q$ curve
!             p(5): high q value of $a, q$ curve
!             p(6): minimum value of q
!     a     : Value of semimajor axis
!
! OUTPUT
!     scat_q : Random value of q (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) seed
!f2py intent(in) np
!f2py intent(in) p
!f2py intent(in) a
    implicit none

    integer (kind=4), intent(inout) :: seed
    integer (kind=4), intent(in) :: np
    real (kind=8), intent(in) :: p(*), a
    integer (kind=4) :: i
    real (kind=8), save :: random, a1, q1, a2, q2, alpha, al, AA, pq, xq, qmax, qmin
    logical, save :: first

    data al /1.3d0/, a1 /43.0d0/, q1 /33.0d0/, a2 /425.0d0/, q2 /46.0d0/, qmin /10.0d0/
    data first /.true./

    if (first) then
       if (np .ge. 1) then
          al = p(1)
       end if
       if (np .ge. 2) then
          a1 = p(2)
       end if
       if (np .ge. 3) then
          q1 = p(3)
       end if
       if (np .ge. 4) then
          a2 = p(4)
       end if
       if (np .ge. 5) then
          q2 = p(5)
       end if
       if (np .ge. 6) then
          qmin = p(6)
       end if
       AA = 10.0**(-al)
       alpha = log10(q1/q2)/log10(a1/a2)
       first = .false.
    end if

    qmax = q1*(a/a1)**alpha
    random=ran_3(seed)
    pq = (1 - AA)*random + AA
    xq = 1.0d0 + log10(pq)/al
    scat_q = xq*(qmax - qmin) + qmin

    return
  end function scat_q

  real (kind=8) function scat_q_2 (seed, np, p, a)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine draws randomly a number according to the scattering q
! distribution, described below.
! The first values have been eye-balled on the scattering q distribution.
!
! The extend of the $q$ distribution depends on $a$. So we first
! scale according to $a$. The minimum value of $q$ is $q_{min} = 10$. The maximum
! value of $q$ is $q_{max} \propto a^\alpha$. From the detectiosn, I determined $a_1
! = 43$ corresponds to $q_1 = 33$ and $a_2 = 425$ corresponds to $q_2 = 46$. This
! sets $\alpha$ (which could be adjusted) and then the proportionnality factor. now
! I define a variable that runs in [0; 1]:
! \begin{displaymath}
! \frac{xq - xq_0}{xq_1 - xq_0}= \frac{q-q_0}{q_1 \left(\frac{a}{a_1}\right)^\alpha - q_0}
! \end{displaymath}
! where $xq_0 = 0$, $xq_1 = 1$, $q_0 = q_{min}$. Now, the distribution of $xq$
! looks like two uniform distributions next to each-other. I'll draw uniformly
! in range [0; xq_b] with fraction fr, and for the rest, uniformly in range
! [xq_b; 1].
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : September 2023
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     np    : Number of parameters (I4)
!     p     : Parameters for distribution (n*R8)
!             p(1): break point $xq_b$ of xq distribution
!             p(2): fraction below xq_b
!             p(3): low a value of $a, q$ curve
!             p(4): low q value of $a, q$ curve
!             p(5): high a value of $a, q$ curve
!             p(6): high q value of $a, q$ curve
!             p(7): minimum value of q
!     a     : Value of semimajor axis
!
! OUTPUT
!     scat_q_2: Random value of q (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) seed
!f2py intent(in) np
!f2py intent(in) p
!f2py intent(in) a
    implicit none

    integer (kind=4), intent(inout) :: seed
    integer (kind=4), intent(in) :: np
    real (kind=8), intent(in) :: p(*), a
    integer (kind=4) :: i
    real (kind=8), save :: random, a1, q1, a2, q2, alpha, xqb, fr, xq, qmax, qmin
    logical, save :: first

    data xqb /0.8d0/, fr /0.4d0/, a1 /43.0d0/, q1 /33.0d0/, a2 /425.0d0/, &
         q2 /46.0d0/, qmin /10.0d0/
    data first /.true./

    if (first) then
       if (np .ge. 1) then
          xqb = p(1)
       end if
       if (np .ge. 2) then
          fr = p(2)
       end if
       if (np .ge. 3) then
          a1 = p(3)
       end if
       if (np .ge. 4) then
          q1 = p(4)
       end if
       if (np .ge. 5) then
          a2 = p(5)
       end if
       if (np .ge. 6) then
          q2 = p(6)
       end if
       if (np .ge. 7) then
          qmin = p(7)
       end if
       alpha = log10(q1/q2)/log10(a1/a2)
       first = .false.
    end if

    qmax = q1*(a/a1)**alpha
    random=ran_3(seed)
    if (random .le. fr) then
       xq = xqb*random/fr
    else
       xq = xqb + (1.0d0-xqb)*(random-fr)/(1.0d0-fr)
    end if
    scat_q_2 = xq*(qmax - qmin) + qmin

    return
  end function scat_q_2

  real (kind=8) function scat_q_3 (seed, np, p, a)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine draws randomly a number according to the scattering q
! distribution, described below.
! The first values have been eye-balled on the scattering q distribution.
!
! The extend of the $q$ distribution depends on $a$. So we first
! scale according to $a$. The minimum value of $q$ is
! $q_{min} = max(10, q_3 (a/a_3)^\beta)$, with $\beta = log10(q_3/q_4)/log10(a_3/a_4)$,
! $a_3 = 100$, $q_3 = 10$, $a_4 = 700$ and $q_4 = 40$
! Given this,
! - $P(q)$ is uniform between 10 and 20 with total fraction fr_1; if q_{min} > 10,
!   then fr_1 is reduced according to reduced range
! - $P(q)$ is uniform between 20 and 32 with total fraction fr_2; if q_{min} > 20,
!   then fr_2 is reduced according to reduced range
! - for the rest, $P(xq)$ is uniform with remaining fraction; xq in range
!   [fr_1+fr_2; 1] corresponds linearly to q in range [max(32, q_{min}); q_{max}]
! with:
! \begin{displaymath}
! \frac{xq - xq_0}{xq_1 - xq_0}= \frac{q-q_0}{q_1 \left(\frac{a}{a_1}\right)^\alpha - q_0}
! \end{displaymath}
! where $xq_0 = fr_1 + fr_2$, $xq_1 = 1$, $q_0 = q_{min}$.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : September 2023
!             In this version I set the various points that define $\alpha$
!             and $\beta$. I can vary fr_1 and fr_2 only
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     np    : Number of parameters (I4)
!     p     : Parameters for distribution (n*R8)
!             p(1): fr_1
!             p(2): fr_2
!     a     : Value of semimajor axis
!
! OUTPUT
!     scat_q_3: Random value of q (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) seed
!f2py intent(in) np
!f2py intent(in) p
!f2py intent(in) a
    implicit none

    integer (kind=4), intent(inout) :: seed
    integer (kind=4), intent(in) :: np
    real (kind=8), intent(in) :: p(*), a
    integer (kind=4) :: i
    real (kind=8), save :: random, a1, qa1, a2, qa2, a3, qa3, a4, qa4
    real (kind=8), save :: q0, q1, q2, qm
    real (kind=8), save :: alpha, beta, fr1, fr2, xq, qmax, qmin, f1, f2
    logical, save :: first

    data a1 /43.0d0/, qa1 /33.0d0/, a2 /425.0d0/, qa2 /46.0d0/, &
         a3 /100.0d0/, qa3 /10.0d0/, a4 /700.0d0/, qa4 /40.0d0/
    data q0 /10.0d0/, q1 /20.0d0/, q2 /32.0d0/
    data fr1 /0.1d0/, fr2 /0.2d0/
    data first /.true./

    if (first) then
       if (np .ge. 1) then
          fr1 = p(1)
       end if
       if (np .ge. 2) then
          fr2 = p(2)
       end if
       alpha = log10(qa1/qa2)/log10(a1/a2)
       beta = log10(qa3/qa4)/log10(a3/a4)
       first = .false.
    end if

    qmax = qa1*(a/a1)**alpha
    qmin = max(q0, qa3*(a/a3)**beta)
    f1 = min(fr1, max(0.0d0, fr1*(q1-qmin)/(q1-q0)))
    f2 = min(fr2, max(0.0d0, fr2*(q2-qmin)/(q2-q1)))
    random=ran_3(seed)
    if (random .le. f1) then
       scat_q_3 = q0 + random*(q1-q0)/f1
    else if (random .le. f1+f2) then
       scat_q_3 = q1 + (random-f1)*(q2-q1)/f2
    else
       qm = max(q2, qmin)
       scat_q_3 = qm + (random-f1-f2)*(qmax-qm)/(1.0d0-f1-f2)
    end if

    return
  end function scat_q_3

  real (kind=8) function scat_q_4 (seed, np, p, a)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine draws randomly a number according to the scattering q
! distribution, described below.
! The first values have been eye-balled on the scattering q distribution.
!
! The extend of the $q$ distribution depends on $a$. So we first
! scale according to $a$. The minimum value of $q$ is
! $q_{min} = max(10, q_3 (a/a_3)^\beta)$, with $\beta = log10(q_3/q_4)/log10(a_3/a_4)$,
! $a_3 = 100$, $q_3 = 10$, $a_4 = 700$ and $q_4 = 40$
! Given this,
! - $P(q)$ is uniform between 10 and 20 with total fraction fr_1; if q_{min} > 10,
!   then fr_1 is reduced according to reduced range
! - $P(q)$ is uniform between 20 and 32 with total fraction fr_2; if q_{min} > 20,
!   then fr_2 is reduced according to reduced range
! - for the rest, $P(xq)$ is increasing as a power of xq with remaining fraction;
!   xq in range [fr_1+fr_2; 1] corresponds linearly to q in range
!   [max(32, q_{min}); q_{max}]
! with:
! \begin{displaymath}
! \frac{xq - xq_0}{xq_1 - xq_0}= \frac{q-q_0}{q_1 \left(\frac{a}{a_1}\right)^\alpha - q_0}
! \end{displaymath}
! where $xq_0 = 0$, $xq_1 = 1$, $q_0 = q_{min}$.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : September 2023
!             In this version I set the various points that define $\alpha$
!             and $\beta$. I can vary fr_1 and fr_2 only
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : Random number generator seed (I4)
!     np    : Number of parameters (I4)
!     p     : Parameters for distribution (n*R8)
!             p(1): fr_1
!             p(2): fr_2
!             p(2): al, the exponent of the increasing distribution of xq
!     a     : Value of semimajor axis
!
! OUTPUT
!     scat_q_4: Random value of q (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in,out) seed
!f2py intent(in) np
!f2py intent(in) p
!f2py intent(in) a
    implicit none

    integer (kind=4), intent(inout) :: seed
    integer (kind=4), intent(in) :: np
    real (kind=8), intent(in) :: p(*), a
    integer (kind=4) :: i
    real (kind=8), save :: random, a1, qa1, a2, qa2, a3, qa3, a4, qa4
    real (kind=8), save :: q0, q1, q2, qm, al, al1
    real (kind=8), save :: alpha, beta, fr1, fr2, xq, qmax, qmin, f1, f2
    logical, save :: first

    data a1 /43.0d0/, qa1 /33.0d0/, a2 /425.0d0/, qa2 /46.0d0/, &
         a3 /100.0d0/, qa3 /10.0d0/, a4 /700.0d0/, qa4 /40.0d0/
    data q0 /10.0d0/, q1 /20.0d0/, q2 /32.0d0/
    data fr1 /0.1d0/, fr2 /0.2d0/
    data first /.true./

    if (first) then
       if (np .ge. 1) then
          fr1 = p(1)
       end if
       if (np .ge. 2) then
          fr2 = p(2)
       end if
       if (np .ge. 3) then
          al = p(3)
       end if
       alpha = log10(qa1/qa2)/log10(a1/a2)
       beta = log10(qa3/qa4)/log10(a3/a4)
       al1 = 1.0d0/al
       first = .false.
    end if

    qmax = qa1*(a/a1)**alpha
    qmin = max(q0, qa3*(a/a3)**beta)
    f1 = min(fr1, max(0.0d0, fr1*(q1-qmin)/(q1-q0)))
    f2 = min(fr2, max(0.0d0, fr2*(q2-qmin)/(q2-q1)))
    random=ran_3(seed)
    if (random .le. f1) then
       scat_q_4 = q0 + random*(q1-q0)/f1
    else if (random .le. f1+f2) then
       scat_q_4 = q1 + (random-f1)*(q2-q1)/f2
    else
       xq = ((random-f1-f2)/(1.0d0-f1-f2))**al1
       qm = max(q2, qmin)
       scat_q_4 = qm + xq*(qmax-qm)
    end if

    return
  end function scat_q_4

  real (kind=8) function qhot (np, p, q)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! This function returns the unnormalized inclination "probability" of having
! perihelion distance "q"
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : April 2020
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! INPUT
!     np    : Number of parameters describing the distribution function (I4)
!     p     : Array of parameters (np*R8)
!     q     : Perihelion distance (R8)
!
! OUTPUT
!     qhot  : Probability of q
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! Set of F2PY directives to create a Python module
!f2py intent(in) np
!f2py intent(in) p
!f2py intent(in) q
    implicit none

! Calling arguments
    integer (kind=4), intent(in) :: np
    real (kind=8), intent(in) :: p(*), q

    qhot = 1./((1.+exp((p(3)-q)/p(4)))*(1.+exp((q-p(1))/p(2))))

    return
  end function qhot

  real (kind=8) function det_q_2(seed, param)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This function draws "q" according 2 uniform distributions over 2
! consecutive ranges, with given fraction in each range.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : April 2023
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : seed for the random number generator (I4)
!     param : parameters for the distribution (4*R8)
!             param(1) = q_min
!             param(2) = q_split
!             param(3) = q_max
!             param(4) = f_low
!
! OUTPUT
!     det_q_2: Random value of q (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! Set of F2PY directives to create a Python module
!f2py intent(in,out) seed
!f2py intent(in) param
    implicit none

    integer (kind=4), intent(inout) :: seed
    real (kind=8), intent(in) :: param(4)
    real (kind=8) :: random

    random = ran_3(seed)
    if (random .le. param(4)) then
! This is the low-q part
       det_q_2 = param(1) + (param(2)-param(1))*random/param(4)
    else
! This is the high-q part
       det_q_2 = param(2) + (param(3)-param(2))*(random-param(4))/(1.0d0-param(4))
    end if

    return
  end function det_q_2

  real (kind=8) function det_q_3(nparam, param, q)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This routine returns the unnormalized q-distribution "probability"
! density for outer/detached objects.
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : April 2023
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     nparam: Number of parameters (I4)
!     param : Parameters (n*R8)
!             param(1) = q_min
!             param(2) = q_split1
!             param(3) = q_split2
!             param(4) = q_max
!     q     : Perihelion distance [au] (R8)
!
! OUPUT
!     det_q_3: Value of the probability (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!f2py intent(in) nparam
!f2py intent(in), depend(nparam) :: param
!f2py intent(in) q
    implicit none

    integer (kind=4), intent(in) :: nparam
    real (kind=8), intent(in) :: param(*), q
    real (kind=8), parameter :: Pi = 3.141592653589793238d0, &
         PiThird = Pi/3.0d0

    if (nparam .ne. 4) stop
    if (q .lt. param(1)) then
       det_q_3 = 0.d0
    else if (q .lt. param(2)) then
       det_q_3 = 0.5d0*(q-param(1))/(param(2)-param(1))
    else if (q .le. param(3)) then
       det_q_3 = 0.5d0
    else if (q .le. param(4)) then
       det_q_3 = 1.d0 - cos((param(4)-q)/(param(4)-param(3))*PiThird)
    else
       det_q_3 = 0.d0
    end if

    return
  end function det_q_3

  real (kind=8) function det_q_n(seed, np, param)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! This function draws "q" according to np-1 uniform distributions over np-1
! consecutive ranges, with given fraction in each range.
!
! BEWARE: sum(param(2:np:2) MUST be 1.0d0
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!
! J-M. Petit  Observatoire de Besancon
! Version 1 : June 2024
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! INPUT
!     seed  : seed for the random number generator (I4)
!     np    : number of parameters    
!     param : parameters for the distribution (4*R8)
!             param(1) = q_min
!             param(2) = fraction of range 1
!             param(3) = end of range 1
!             param(4) = fraction of range 2
!             param(5) = end of range 2
!             ...
!
! OUTPUT
!     det_q_2: Random value of q (R8)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
! Set of F2PY directives to create a Python module
!f2py intent(in,out) seed
!f2py intent(in) np
!f2py intent(in) param
    implicit none

    integer (kind=4), intent(inout) :: seed
    integer (kind=4), intent(in) :: np
    real (kind=8), intent(in) :: param(np)
    real (kind=8) :: random, s1
    integer (kind=4) :: i

    random = ran_3(seed)
    s1 = 0.0d0
    do i = 2, np, 2
       if ((random .ge. s1) .and. (random .le. s1+param(i))) then
          det_q_n = param(i-1) + (param(i+1)-param(i-1))*(random-s1)/param(i)
          exit
       end if
       s1 = s1 + param(i)
    end do

    return
  end function det_q_n

  real (kind=8) function ran_3(idum)
!f2py intent(in,out) idum
    INTEGER (KIND=4), intent(inout) :: idum
    INTEGER (KIND=4), parameter :: MBIG=1000000000, MSEED=161803398, MZ=0
    REAL (KIND=8), parameter :: FAC=1.d0/MBIG
    INTEGER :: i,ii,k,mj,mk
    INTEGER (KIND=4), save :: iff,inext,inextp,ma(55)
    data iff /0/
    if(idum.lt.0.or.iff.eq.0)then
       iff=1
       mj=abs(MSEED-abs(idum))
       mj=mod(mj,MBIG)
       ma(55)=mj
       mk=1
       do i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
       end do
       do k=1,4
          do i=1,55
             ma(i)=ma(i)-ma(1+mod(i+30,55))
             if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
          end do
       end do
       inext=0
       inextp=31
       idum=1
    endif
    inext=inext+1
    if(inext.eq.56)inext=1
    inextp=inextp+1
    if(inextp.eq.56)inextp=1
    mj=ma(inext)-ma(inextp)
    if(mj.lt.MZ)mj=mj+MBIG
    ma(inext)=mj
    ran_3=mj*FAC
    return
  end function ran_3

end module modelutils
