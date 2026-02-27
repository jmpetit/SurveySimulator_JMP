module gimeobjut

  use datadec
  use elemutils
  use rot
  use modelutils
  use ioutils
  use xvutils
  use numutils

!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! First define variables so they are accessible from a Python wrapper
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! Values of time and planet positions for all models
  real (kind=8) :: lambdaN, epoch_m
  common /com_time/ epoch_m, lambdaN
  data lambdaN /5.489d0/, epoch_m /2453157.5d0/

contains
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! Generic model routines for Survey Simulator, version 2.0 for OSSOS
!
! Calling sequence in SurveySimulator.f (survey simulator driver) is:
!
! Loop (some condition):
!     call GiMeObj(arg_list_1)
!     Check model ended:
!         set exit condition
!     call Detos1(arg_list_2)
!     Check detection and tracking:
!         store results
!
! where arg_list_1 is
! (filena, seed, a, e, inc, node, peri, M, epoch, h, color, gb, ph,
!  period, amp, commen, nchar, ierr)
! with:
!
! INPUT
!     filena: name of containing description of model, read in by the
!             model subroutine "modname" (CH)
!     seed  : Random number generator seed (I4)
!
! OUTPUT
!     a     : semimajor axis (R8)
!     e     : eccentricity (R8)
!     inc   : Inclination [rad] (R8)
!     node  : Longitude of node [rad] (R8)
!     peri  : Argument of perihelion [rad] (R8)
!     M     : Mean anomaly [rad] (R8)
!     epoch : epoch of the orbital elements, in Julian Day (R8)
!     h     : absolute magnitude of object in band filter "x" (R8)
!     color : array of colors "y-x", where the index of "y" is as
!             described in detos1 (10*R8)
!                color(1) : g-x
!                color(2) : r-x
!                color(3) : i-x
!                color(4) : z-x
!                color(5) : u-x
!                color(6) : V-x
!                color(7) : B-x
!                color(8) : R-x
!                color(9) : I-x
!     gb    : opposition surge factor, Bowell formalism (R8)
!     ph    : phase of lightcurve at epoch [rad] (R8)
!     period: period of lightcurve [day] (R8)
!     amp   : amplitude of lightcurve [mag] (R8)
!     commen: user specified string containing whatever the user wants (CH*100)
!     nchar : number of characters in the comment string that should be
!             printed out in output files if the object is detected;
!             maximum of 100 (I4)
!     ierr  : return code
!                  0 : nominal run, things are good
!                100 : end of model, exit after checking this object
!                -10 : could not get all orbital elements, skip object
!                -20 : something went grossly wrong, should quit
!
! The model subroutines can access files using logical unit numbers from
! 10 to 15. This range in reseved for them and won't be used by the
! drivers nor SurveySubs routines.
!
! It is good practice that when first started, the GiMeObj routine
! writes a file describing the model used, the versino and the date of
! the routine.
!
! Since this routine is called once for every object created, it needs
! to get all the required parameters once when it is called the first
! time, then save these values for future use.
!
! The following routine gives a working example of a model routine. It
! is probably worth reading it through.
!
! The survey simulator expects orbital elements with respect to ecliptic
! reference frame.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! File generated on 2023-07-31
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  subroutine GiMeObj (filena, seed, o_m, epoch, h, color, gb, ph, period, &
       amp, commen, nchar, ierr)
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! This routine generates an object from a model parametric model of the
! inner population.
!
! Version 1.0 draws (a, q) according to CFEPS model, i according to =hot_inc=
! distribution with respect to forced reference frame, H distribution is
! =H_dist_hot_3= with possible knee or divot at some mag.
!
! Version 1.1 draws (a, q) according to CFEPS model, i according to Brown
! distribution with respect to forced reference frame, H distribution is
! =H_dist_hot_3= with possible knee or divot at some mag.
!
! Version 1.2 draws (a, q) according to CFEPS model, i according to Brown
! distribution with respect to forced reference frame, H distribution is
! =H_dist_hot_4= with possible knee or divot at some mag.
!
! Version 1.3 draws (a, q) according to CFEPS model, i according to Brown
! distribution with respect to forced reference frame, H distribution is
! =H_draw_hot_6= with possible knee or divot at some mag.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! J-M. Petit  Observatoire de Besancon
! Version 1.0 : July 2023
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! INPUT
!     filena: name of containing description of model, read in by the
!             model subroutine "modname" (CH)
!     seed  : Random number generator seed (I4)
!
! OUTPUT
!     o_m   : orbital elements of object (t_orb_m)
!     epoch : Time of elements [JD] (R8)
!     h     : Absolute magnitude of object in 'x' band, what ever this is (R8)
!     color : Array of colors (10*R8)
!                colors(1) : g-x
!                colors(2) : r-x
!                colors(3) : i-x
!                colors(4) : z-x
!                colors(5) : u-x
!                colors(6) : V-x
!                colors(7) : B-x
!                colors(8) : R-x
!                colors(9) : I-x
!     gb    : opposition surge factor, Bowell formalism (R8)
!     ph    : phase of lightcurve at epoch [rad] (R8)
!     period: period of lightcurve [day] (R8)
!     amp   : amplitude of lightcurve [mag] (R8)
!     commen: user specified string containing whatever the user wants (CH*100)
!     nchar : number of characters in the comment string that should be
!             printed out in output files if the object is detected;
!             maximum of 100 (I4)
!     ierr  : return code (I4)
!                  0 : nominal run, things are good
!                100 : end of model, exit after checking this object
!                -10 : could not get all orbital elements, skip object
!                -20 : something went grossly wrong, should quit
!
! The user can fill the 100-character 'commen' string any way they
! wish; this comment string will be printed by the driver on the output
! line of each detection.  Examples of the comment might be resonance name
! and libration amplitude, or the name of a component in the GiMeObj model
! that the object responds to. The nchar variable (passed back to Driver)
! allows the user to print only the first nchar characters of this string.
!
! This routine uses logical unit 10 to access the file containing the model.
!
! The model uses the following indices in incdism to define the various
! distributions. Remember that only indices from 1 to 10 are allowed.
!    1:
!    2:
!    3:
!    4:
!    8:
!    9:
!   10:
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!f2py intent(in) filena
!f2py intent(in,out) seed
!f2py intent(out) o_m
!f2py intent(out) epoch
!f2py intent(out) h
!f2py intent(out) color
!f2py intent(out) gb
!f2py intent(out) ph
!f2py intent(out) period
!f2py intent(out) amp
!f2py intent(out) commen
!f2py intent(out) nchar
!f2py intent(out) ierr
    implicit none

! Calling arguments
    integer (kind=4), intent(inout) :: seed
    integer (kind=4), intent(out) :: ierr, nchar
    type(t_orb_m), intent(out) :: o_m
    real (kind=8), intent(out) :: epoch, h, color(10), gb, ph, period, amp
    character(*), intent(in) :: filena
    character(100), intent(out) :: commen

! Some values better set up as parameters
    integer, parameter :: &
         lun_m = 20,            &! Logical unit number for data file reading
         lun_ll = 21             ! Logical unit number for logging
    real (kind=8), parameter :: &
         Pi = 3.141592653589793238d0, &! Pi
         TwoPi = 2.0d0*Pi,      &! 2*Pi
         drad = Pi/180.0d0       ! Degree to radian convertion: Pi/180

! Internal storage
    type(t_v3d) :: p, opos, ovel
    type(func_holder) :: fp
    character(5) :: zone         ! Time zone
    character(8) :: date         ! Date of execution
    character(10) :: time        ! Time of execution
    integer (kind=4), save :: &
         values(8),             &! Date and time of execution
         flag,                  &! Tell invar_ecl_osc which direction to go
         i,                     &! Dummy index
         comp,                  &! Index of the component
         n_h,                   &! Number of parameters for H distributions
         nlog,                  &! Number of objects to log
         nlogged,               &! Number of objects already logged
         nparam                  ! Number of parameters for the function called
                                 ! by routine incdism
    real (kind=8), save :: &
         hcut,                  &! Largest Hx value in model
         h_params(60),          &! Parameters for the H distribution
         inc_f, node_f, peri_f, &! Free inclination and node and arg of peri
         i_ref, om_ref,         &! Coordinates of forced plane
         param(50),             &! Temporary storage for distribution parameters
         q,                     &! Perihelion distance
         random,                &! Random number
         r,                     &! Distance of object to Sun
         colorH(10),            &! Color parameters of model for hot component
         ra, dec,               &!
         delta, or, mag, alpha, &!
         inimin,                &! Lower limit of secular instability zone
         inimax,                &! Upper limit of secular instability zone
         beta_ah,               &! Index of the a distribution
         ah_min, ah_max,        &! Lower and upper limit of a distribution
         a0sl, a1sl              ! Intermediate values for a distribution
    logical, save :: &
         first                   ! Tells if first call to routine
    real (kind=8) :: &
         epoch_m,               &! Epoch of elements [JD]
         lambdaN                 ! Longitude of Neptune at epoch
! Lightcurve and opposition surge effect parameters
    real (kind=8), save :: gb0, ph0, period0, amp0

! Place some variables in common block so they can be accessed directly
! by a Python program.
    common /com_time/ epoch_m, lambdaN

! Sets initial values
    data &
         first /.true./,        &! First call
         gb0     / 0.15d0/,     &! Opposition surge effect
         ph0     / 0.00d0/,     &! Initial phase of lightcurve
         period0 / 0.60d0/,     &! Period of lightcurve
         amp0    / 0.00d0/       ! Amplitude of lightcurve (peak-to-peak)

! This is the first call
    if (first) then
! Define the region of the inner belt
       inimin = 7.d0  *drad   ! This cuts out 7-20 deg inclinations (nu 8)
       inimax = 20.d0 *drad
! Reads in other parameters describing the model.
       open (unit=lun_m, file=filena, status='old', err=1000)
       read (lun_m, *) n_h
! hot
       read (lun_m, *) (h_params(i+0*n_h), i=1,n_h)
       comp = 1
       read (lun_m, *) (param(comp*10+i),i=1,2)
       do i = 1, 2
          param(comp*10+i) = param(comp*10+i)*drad
       end do
       comp = 2
       read (lun_m, *) (param(comp*10+i),i=1,4)
       read (lun_m, *) beta_ah, ah_min, ah_max
       a0sl = ah_min**(beta_ah + 1.d0)
       a1sl = ah_max**(beta_ah + 1.d0)
       read (lun_m, *) (colorH(i),i=1,10)      ! read color array
       read (lun_m, *) nlog
       close (lun_m)
! Writes a file describing the model that was used.
       open (unit=lun_ll, file='ModelUsed.dat', access='sequential', &
            status='unknown')
       write (lun_ll, '(a)') '# Inner population model.'
       write (lun_ll, '(a)') '# Version OSSOS 8.1, 2023-07-31'
       call date_and_time(date, time, zone, values)
       write (lun_ll, '(a17,a23,2x,a5)') '# Creation time: ', &
            date(1:4)//'-'//date(5:6)//'-'//date(7:8)//'T' &
            //time(1:2)//':'//time(3:4)//':'//time(5:10), zone
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a,10(1x,f5.2))') &
            '# Colors for hot: ', (colorH(i),i=1,10)
       write (lun_ll, '(''#'')')
       comp = 1
       write (lun_ll, '(a,2(1x,f5.2))') &
            '# Parameters for inclination distribution: ', &
            (param(comp*10+i)/drad, i=1,2)
       comp = 2
       write (lun_ll, '(a,4(1x,f5.2))') &
            '# Parameters for q distribution: ', &
            (param(comp*10+i), i=1,4)
       write (lun_ll, '(a,3(1x,f6.2))') &
            '# Parameters for a distribution: ', &
            beta_ah, ah_min, ah_max
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a,4(1x,f5.2))') &
            '# Hot H-dist. parameter:          ', (h_params(i+0*n_h), i=1,n_h)
       write (lun_ll, '(''#'')')
       write (lun_ll, '(a,a,a,a)') &
            '#   a        e        i      Omega    omega      M', &
            '        H      dist   comment     ifree   Omfree ', &
            '     epoch      omfree     iref   Omref   ', &
            '   RA        DEC    alpha  mag'
       close (lun_ll)
       call ObsPos(500, epoch_m, opos, ovel, or, ierr)

       nlogged = 0

! Change "first" so this is not called anymore
       first = .false.
    end if
!
! Only 1 component, inner
    commen = 'inner_'
    nchar = 9

!
! Draw inclination
1150 continue
!    nparam = 4
!    fp%f_ptr => hot_inc
!    call incdism (seed, nparam, param(10+1), 0.0d0*drad, 50.d0*drad, inc_f, &
!         1, ierr, fp)
    nparam = 2
    fp%f_ptr => offgau
    call incdism (seed, nparam, param(10+1), 0.0d0*drad, 70.d0*drad, inc_f, &
         3, ierr, fp)
! Rejects if inclination in the secular instability zone
    if ((inc_f .gt. inimin) .and. (inc_f .lt. inimax)) goto 1150

1100 continue
!
! Draw 'a'
    random=ran_3(seed)
    o_m%a = (a0sl + (a1sl-a0sl)*random)**(1.0d0/(beta_ah+1.0d0))
!
! Index of distributions for incdism
!   1: hot_inc
!   2: qhot
!   3: offgau
!   4: h_dist_hot_4
!
! Draw 'q'
    nparam = 4
    fp%f_ptr => qhot
    call incdism (seed, nparam, param(20+1), 33.0d0, ah_max, q, &
         2, ierr, fp)
! Physical limitaion on q
    if (q .ge. o_m%a) goto 1100
! Put in a cut for instability at low q and low i
    if (q .lt. 37.d0-inc_f/drad*0.2d0) goto 1150
!
! Determination of "e"
    o_m%e = 1.d0 - q/o_m%a
!
! H-mag distribution: Exponential cutoff and broken exponential law
! Also set colors for object
    h = H_draw_hot_6(seed, n_h, h_params(0*n_h+1))
    do i = 1, 10
       color(i) = colorH(i)
    end do
!
! Angles: uniform distribution on allowable values
    random=ran_3(seed)
    o_m%node = random*TwoPi
    random=ran_3(seed)
    o_m%peri = random*TwoPi
    random=ran_3(seed)
    o_m%m = random*TwoPi
!
! Set up epoch for orbial elements
    epoch = epoch_m

! Define values for lightcurve and opposition surge effect
    gb = gb0
    ph = ph0
    period = period0
    amp = amp0
!
! The model above gives orbital elements with respect to the forced
! plane reference frame (orientation depending on 'a')
! The survey simulator expects the orbital elements with respect to the
! ecliptic, so convert them.
    o_m%inc = inc_f
    node_f = o_m%node
    peri_f = o_m%peri
    call forced_plane_damp(o_m%a, inc_f/drad, i_ref, om_ref)
    flag = 1
    call ref_ecl_osc (flag, o_m, o_m, i_ref*drad, om_ref*drad, ierr)
!
! Store object if user requested
! Normally, we should explicitly open a file and write to its end, it
! seems like the pointer to the file is not retained from one call to
! the other, so simply use the default file assigned to the logical
! unit. In this case, the output file will be something like "fort.11"
    if ((nlog .lt. 0) .or. (nlogged .lt. nlog)) then
       call pos_cart(o_m, p)
       call DistSunEcl(epoch, p, r)
       call RADECeclXV(p, opos, delta, ra, dec)
       call AppMag(r, delta, or, h, gb, alpha, mag, ierr)
       open (unit=lun_ll, file='ModelUsed.dat', access='append', &
            status='old')
       write(lun_ll,101) o_m%a, o_m%e, o_m%inc/drad, o_m%node/drad, &
            o_m%peri/drad, o_m%m/drad, h, sqrt(r*delta), commen(1:nchar), &
            inc_f/drad, node_f/drad, epoch, peri_f/drad, &
            i_ref, om_ref, ra/drad, dec/drad, alpha/drad, mag
       close (lun_ll)
101    format(f9.4,1x,5(f8.4,1x),f6.2,1x,f9.4,1x,a9,2(1x,f8.4), &
            1x,f13.5,3(1x,f8.4),1x,f9.5,1x,f9.5,1x,f5.2,1x,f5.2)
       nlogged = nlogged + 1
    end if

! Prepare return code
    ierr = 0

    return

1000 continue
! If we get here, there is something really wrong, better return with
! panic code.
    ierr = -20
    return

  end subroutine GiMeObj

end module gimeobjut
