program DriverPop
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
! File DriverPop.f95
!
! J.-M. Petit  Observatoire de Besac
! Version 1:  May 2013
! Version 2:  April 2023
!
! The purpose of this program is to compare models of the Kuiper Belt to
! reality by actually comparing what the surveys (say CFEPS or OSSOS)
! would have found if the model was a good representation of the real
! world, to the objects that were actually found.
!
! This version of the Driver will also write a file to compute the calibrated
! population at a given H mag by running 'n_pop' times the detection of 'n_real'
! simulated detections
!
! The workflow of the survey simulator driver is as follows:
!
!     Loop (some condition):
!         call GiMeObj(arg_list_1)
!         Check model ended:
!             set exit condition
!         call Detos1(arg_list_2)
!         Check detection and tracking:
!             store results
!
! The GiMeObj routine is in charge of providing a single new object at
! each call. The Detos1 routine determines is the proposed object would
! have been detected by the survey.
!
! The survey simulator expects orbital elements with respect to ecliptic
! reference frame.
!
! Logical unit numbers 7 to 19 are reserved to use by Driver.f and
! SurveySubs.f and should not be used by GiMeObj or any other routine
! that you may add to the driver. Please use logical unit numbers
! starting from 20.
!
! As currently written, the driver module 'gimeobjut'
! containing the definition of the model. The correct way to use this
! feature is to have one's GiMeObj routine in a file <whatever.f> and
! create a symobolic link: 
!
!     ln -s <whatever.f> GiMeObj.f
!
! For more information, please refer to README files in the source directory.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

  use gimeobjut
  use surveysub

  implicit none

  integer, parameter :: n_obj_max = 10000, screen = 6, keybd = 5, verbose = 9
  integer :: lun_h, lun_t, lun_p
  type(t_orb_m) :: o_m
! color array NEEDS to be length 10 or more!
  real (kind=8) :: h, epoch, m_int, d_ra, d_dec, r, delta, ra, dec, random, &
       mt, color(10), gb, ph, period, amp, jday_p, m_rand, eff, rn_iter, &
       rn0, eff_lim, h_rand, h_pop, rn_i_h
  integer :: n_hits, n_track, n_h_h, n_t_h, n_i_h, ierr, seed, flag, isur, ic, &
       n_iter, n_track_max, nchar, values(8), i1, i2, nh0, nt0, n_real, n_pop
  character(80) :: distri_file, trk_outfile, det_outfile, pop_outfile, line
  character(100) :: survey_dir, comments
  character(10) :: surna, time
  character(8) :: date
  character(5) :: zone
  logical :: keep_going, finished

  lun_h = 10
  lun_t = 11
  lun_p = 12
  keep_going = .true.

! Get arguments
! Seed for random number generator
  read (5, *, err=9999) seed
! Number of real detections
  read (5, *, err=9999) n_real
! Number of times we draw 'n_real' detections
  read (5, *, err=9999) n_pop
  n_real = abs(n_real)
  n_track_max = n_pop*n_real
! H value at which we want the population estimate
  read (5, *) h_pop
! Directory containing the characterization files
  read (5, '(a)', err=9999) survey_dir
  survey_dir = strip_comment(survey_dir)
! File with model parameters
  read (5, '(a)', err=9999) distri_file
  distri_file = strip_comment(distri_file)
! Name for the detected objects outfile
  read (5, '(a)', err=9999) det_outfile
  det_outfile = strip_comment(det_outfile)
! Name for the tracked objects outfile
  read (5, '(a)', err=9999) trk_outfile
  trk_outfile = strip_comment(trk_outfile)
! Name for the population outfile
  read (5, '(a)', err=9999) pop_outfile
  pop_outfile = strip_comment(pop_outfile)

!  print *, n_real, n_pop, n_track_max, h_pop

! Open output files and write header
  open (unit=lun_h, file=det_outfile, status='new', err=9500)
  open (unit=lun_t, file=trk_outfile, status='new', err=9501)
  open (unit=lun_p, file=pop_outfile, status='new', err=9504)

! copy header of model file into trk and detect files
  open (unit=90, file=distri_file, status='old', err=9503)
222 continue
  read (90, '(a)', err=9999) line
  if ( line(1:1) .eq. '#' ) then
     !    don't write out if this is the header line.. hack.
     if (index(line, 'peri') .eq. 0) then 
        write(lun_t, '(a)') line
        write(lun_h, '(a)') line
        write(lun_p, '(a)') line
     endif
     goto 222
  endif
  close (90)

  write (lun_h, '(''# Seed: '', i10)') seed
  write (lun_h, '(''#'')')
  call date_and_time(date, time, zone, values)
  write (lun_h, '(a17,a23,2x,a5)') '# Creation time: ', &
       date(1:4)//'-'//date(5:6)//'-'//date(7:8)//'T'  &
       //time(1:2)//':'//time(3:4)//':'//time(5:10), zone
  write (lun_h, '(''#'')')
  write (lun_h, ' &
       (''# flag: >0: detected; >2: characterized; 0 mod(2): tracked'')')
  write (lun_h, '(''# Survey: name of the block'')')
  write (lun_h, '(''#'')')
  write (lun_h, '(a,a,a,a)') &
       '#   a      e        i       node     peri      M       H_int ', &
       '    q         r      delta', &
       '    mrand    color flag m_int    H_rand    eff   RA(H)  ', &
       '   DEC    Surv. Comments'

  write (lun_t, '(a17,a23,2x,a5)') '# Creation time: ', &
       date(1:4)//'-'//date(5:6)//'-'//date(7:8)//'T'  &
       //time(1:2)//':'//time(3:4)//':'//time(5:10), zone
  write (lun_t, '(a,a)') &
       '#   a           e        i      Node     peri       M       H_int     ', &
       '   q      r       delta  m_rand h_rand color Comments'

  write (lun_p, '(a17,a23,2x,a5)') '# Creation time: ', &
       date(1:4)//'-'//date(5:6)//'-'//date(7:8)//'T'  &
       //time(1:2)//':'//time(3:4)//':'//time(5:10), zone
  write (lun_p, '(a31,f5.2)') '# Population estimate to H_r = ', h_pop
  write (lun_p, '(a)') '#  n_pop      n_det     n_trk'

! Initialize counters
  n_hits = 0
  n_track = 0
  n_iter = 0
  rn_iter = 0.0d0
  n_h_h = 0
  n_t_h = 0
  n_i_h = 0
  rn_i_h = 0.0d0
  rn0 = 0.d0
  nh0 = 0
  nt0 = 0

! Open and read in object distribution
100 continue
      
!     Main loop: loop on objects
  if (keep_going) then

!        Select object.
     nchar = 0
     call GiMeObj (distri_file, seed, o_m, epoch, h, color, gb, ph, period, &
          amp, comments, nchar, ierr)

!     print *, o_m, epoch

     if (ierr .eq. -10) then       ! Something wrong with this object, go to next one
        goto 100
     else if (ierr .eq. -20) then     ! Something very wrong happend, stop
        write (screen, '(a)') 'GiMeObj returned -20, stopping.'
        goto 2010
     else if (ierr .eq. 100) then     ! Reached end of model, prepare to stop
        keep_going = .false.
     end if

!        Determine if the object would be detected
     call Detos1 (o_m, epoch, h, color, gb, ph, period, amp, survey_dir, seed, &
          flag, ra, dec, d_ra, d_dec, r, delta, m_int, m_rand, eff, isur, mt, &
          jday_p, ic, surna, h_rand, ierr)

     if (ierr < 0) then
         write(screen, *) "Failed while attempting to determine if object is detectable: ", ierr
         goto 2011
     end if

!        Count number of iterations; may exceed 2**31, so use a real*8 helper
     n_iter = n_iter + 1
     if (n_iter .gt. 2000000000) then
        rn_iter = rn_iter + dble(n_iter)
        n_iter = 0
        if (rn_iter .ge. 9.d9) goto 2000
     end if
! Count number of iterations with H <= H_pop
     if (h .le. h_pop) n_i_h = n_i_h + 1
     if (n_i_h .gt. 2000000000) then
        rn_i_h = rn_i_h + dble(n_i_h)
        n_i_h = 0
        if (rn_i_h .ge. 9.d9) goto 2000
     end if

!        Check if detected and tracked, and write to output files
     if (flag .gt. 0) then

!        m_int and h are in "x" band (filter of object creation)
!        m_rand and h_rand are in discovery filter
        n_hits = n_hits + 1
! Number of detections that are bright enough
        if (h .le. h_pop) n_h_h = n_h_h + 1
        write (lun_h, 9000) o_m%a, o_m%e, o_m%inc/drad, o_m%node/drad, &
             o_m%peri/drad, o_m%m/drad, h, &
             o_m%a*(1.d0-o_m%e), r, delta, &
             m_rand, color(ic), &
             flag, m_int, h_rand, eff, ra/drad/15., dec/drad, &
             surna, comments(1:nchar)
        if ((flag .gt. 2) .and. (mod(flag,2) .eq. 0)) then
           n_track = n_track + 1
           write (lun_t, 9010) o_m%a, o_m%e, o_m%inc/drad, o_m%node/drad, &
                o_m%peri/drad, o_m%m/drad, h, &
                o_m%a*(1.d0-o_m%e), r, delta, &
                m_rand, h_rand, color(ic), comments(1:nchar)
! Number of tracked objects that are bright enough
           if (h .le. h_pop) then
              n_t_h = n_t_h + 1
              if (mod(n_t_h, n_real) .eq. 0) then
                 write (lun_p, *) rn_i_h + dble(n_i_h) - rn0, &
                      n_h_h - nh0, n_t_h - nt0
                 rn0 = rn_i_h + dble(n_i_h)
                 nh0 = n_h_h
                 nt0 = n_t_h
              end if
           end if
        end if
     end if

! Should we continue ?
     if (n_t_h .ge. n_track_max) then
        keep_going = .false.
     end if
     goto 100
!     end of the if( keep_going ) loop    
  end if

2000 continue
  write (lun_h, '(''#'')')
  write (lun_h, '(''# Total number of objects:   '', f11.0)') &
       rn_iter + dble(n_iter)
  write (lun_h, '(''# Number of detections:      '', i7)') n_hits
  write (lun_h, '(''# Number of tracked objects: '', i7)') n_track
  write (lun_h, '(''#'')')
  write (lun_h, '(''# For objects brighter than Hr = '',f5.2)') h_pop
  write (lun_h, '(''# Total number of objects:   '', f11.0)') &
       rn_i_h + dble(n_i_h)
  write (lun_h, '(''# Number of detections:      '', i7)') n_h_h
  write (lun_h, '(''# Number of tracked objects: '', i7)') n_t_h
  close (lun_h)
  close (lun_t)
  close (lun_p)

  call exit (0)

9000 format (f8.3,1x,f6.3,1x,4(f8.3,1x),f8.2,1x,5(f8.3,1x),i2,2(1x,f8.3), &
          1x,f6.2,1x,f8.5,1x,f8.4,1x,a6,1x,a)
9010 format (6(f9.4,1x),f5.2,1x,3(f9.4,1x),3(f5.2,1x),a)

9503 continue
  write (screen, *) 'File === "', distri_file, '" Does not exit?'
  call exit(-1)

9500 continue
  write (screen, *) 'File --- "', det_outfile, '" already exists. '
  goto 9502

9501 continue
  write (screen, *) 'File +++ "', trk_outfile, '" already exists. '
  goto 9502

9502 continue
  write (screen, *) 'Make sure "', det_outfile, '", "', trk_outfile, &
       '" and "', pop_outfile, '" do not exist and restart SurveySimulator.'
  call exit(-1)

9504 continue
  write (screen, *) 'File %%% "', pop_outfile, '" already exists. '
  goto 9502

9999 continue
  write (screen, *) 'Usage: SurveySimulator < input'
  write (screen, *)
  write (screen, *) 'This will read the following from the keyboard:'
  write (screen, *) '<seed>'
  write (screen, *) '<n_track_max>'
  write (screen, *) '<survey_dir>'
  write (screen, *) '<distrib_file>'
  write (screen, *) '<det_file>'
  write (screen, *) '<track_file>'
  write (screen, *) 'where:'
  write (screen, *) &
       '<seed>: integer used as seed for the random number generator'
  write (screen, *) '<n_track_max>: maximum number of simulated tracked' &
       //' detection if > 0,'
  write (screen, *) '              -maximum number of trials if < 0'
  write (screen, *) '<survey_dir>: directory of survey charactarization files'
  write (screen, *) '<distrib_file>: name of the input file for GiMeObj'
  write (screen, *) '<det_file>: output file for detected objects'
  write (screen, *) '<track_file>: output file for tracked objects'

  call exit(0)

2010 continue
  write (lun_h,*) "# Failed while loading objects from "//distri_file
  call exit(-20)

2011 continue
  write (lun_h,*) "# Failed while loading survey information "
  call exit(-20)


end program DriverPop
