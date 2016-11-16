!----------------------------------------------------------------------------
! This file is part of UCLALES.
!testi
! UCLALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! UCLALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 199-2007, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
module step

  implicit none

  integer :: istpfl = 1
  real    :: timmax = 18000.
  logical :: corflg = .false.
  logical :: rylflg = .true.

  real    :: frqhis =  9000.
  real    :: frqanl =  3600.
  real    :: radfrq =  0.

  real    :: time   =  0.
  real    :: strtim =  0.0
  real    :: cntlat =  31.5 ! 30.0
  logical :: outflg = .true.

contains
  !
  ! ----------------------------------------------------------------------
  ! Subroutine model:  This is the main driver for the model's time
  ! integration.  It calls the routine tstep, which steps through the
  ! physical processes active on a time-step and updates variables.  It
  ! then checks to see whether or not different output options are
  ! satisfied.
  subroutine stepper

    use mpi_interface, only : myid, double_scalar_par_max

    use grid, only : dtl, dzt, zt, zm, nzp, dn0, u0, v0, a_up, a_vp, a_wp, &
         a_uc, a_vc, a_wc, write_hist, write_anal, close_anal, dtlt,  &
         dtlv, dtlong, nzp, nyp, nxp, level,                          &
         ! For mass budged
         a_rp, a_rc, a_srp, a_dn

    use stat, only : sflg, savg_intvl, ssam_intvl, write_ps, close_stat, mcflg, acc_massbudged,  &
         write_massbudged
    use thrm, only : thermo

    logical, parameter :: StopOnCFLViolation = .False.
    real, parameter :: cfl_upper = 0.50, cfl_lower = 0.30

    real    :: t1,t2,tplsdt,begtime,cflmax,gcflmax
    integer :: istp, iret
    logical :: cflflg
    !
    ! Timestep loop for program
    !
    begtime = time
    istp = 0

    do while (time + 0.1*dtl < timmax)

       call cpu_time(t1)           !t1=timing()

       istp = istp+1
       tplsdt = time + dtl + 0.1*dtl
       sflg = (min(mod(tplsdt,ssam_intvl),mod(tplsdt,savg_intvl)) < dtl  &
            .or. tplsdt >= timmax  .or. tplsdt < 2.*dtl)

       call t_step(cflflg,cflmax)

       time  = time + dtl

       call double_scalar_par_max(cflmax,gcflmax)
       cflmax = gcflmax

       if (cflmax > cfl_upper .or. cflmax < cfl_lower) then
          call tstep_reset(nzp,nxp,nyp,a_up,a_vp,a_wp,a_uc,a_vc,a_wc,     &
               dtl,dtlong,cflmax,cfl_upper,cfl_lower)
          dtlv=2.*dtl
          dtlt=dtl
       end if

       !
       ! output control
       !
       if (mod(tplsdt,savg_intvl)<dtl .or. time>=timmax .or. time==dtl)   &
       call write_ps(nzp,dn0,u0,v0,zm,zt,time)

       if ((mod(tplsdt,frqhis) < dtl .or. time >= timmax) .and. outflg)   &
            call write_hist(2, time)
       if (mod(tplsdt,savg_intvl)<dtl .or. time>=timmax .or. time==dtl)   &
            call write_hist(1, time)

       if ((mod(tplsdt,frqanl) < dtl .or. time >= timmax) .and. outflg) then
          call thermo(level)
          call write_anal(time)
       end if

       if (cflflg) then
          cflflg=.False.
          if (StopOnCFLViolation) call write_hist(-1,time)
       end if

       if(myid == 0) then
          call cpu_time(t2)           !t1=timing()
          if (mod(istp,istpfl) == 0 ) print "('   Timestep # ',i5," //     &
              "'   Model time(sec)=',f10.2,3x,'CPU time(sec)=',f8.3)",     &
              istp, time, t2-t1
       endif

    enddo

    IF (mcflg) THEN
       !
       ! Juha:
       ! Get the final statistics of atmospheric water for mass budged
       CALL acc_massbudged(nzp,nxp,nyp,1,dtlt,dzt,a_dn,    &
            rv=a_rp,rc=a_rc,prc=a_srp)

       CALL write_massbudged

       !
       !
    END IF ! mcflg

    call write_hist(1, time)
    iret = close_anal()
    iret = close_stat()

  end subroutine stepper
  !
  !----------------------------------------------------------------------
  ! subroutine tstep_reset: Called to adjust current velocity and reset
  ! timestep based on cfl limits
  !
  subroutine tstep_reset(n1,n2,n3,up,vp,wp,uc,vc,wc,dtl,dtmx,cfl,c1,c2)

  integer, intent (in) :: n1,n2,n3
  real, intent (in)    :: up(n1,n2,n3),vp(n1,n2,n3),wp(n1,n2,n3),dtmx,cfl,c1,c2
  real, intent (inout) :: uc(n1,n2,n3),vc(n1,n2,n3),wc(n1,n2,n3),dtl

  integer :: i,j,k
  real    :: cbar, dtl_old

  cbar = (c1+c2)*0.5
  dtl_old = dtl

  if (cfl > c1) dtl = min(dtmx,dtl*cbar/c1)
  if (cfl < c2) dtl = min(dtmx,dtl*cbar/c2)

  do j=1,n3
     do i=1,n2
        do k=1,n1
           uc(k,i,j) = up(k,i,j) + (uc(k,i,j)-up(k,i,j))*dtl/dtl_old
           vc(k,i,j) = vp(k,i,j) + (vc(k,i,j)-vp(k,i,j))*dtl/dtl_old
           wc(k,i,j) = wp(k,i,j) + (wc(k,i,j)-wp(k,i,j))*dtl/dtl_old
        end do
     end do
  end do

end subroutine tstep_reset

  !
  !----------------------------------------------------------------------
  ! subroutine t_step: Called by driver to timestep through the LES
  ! routines.  Within many subroutines, data is accumulated during
  ! the course of a timestep for the purposes of statistical analysis.
  !
  subroutine t_step(cflflg,cflmax)

    use grid, only : level, dtl, dtlt, Tspinup,                                         &
                     ! Added parameters for interfacing with SALSA
                     nxp, nyp, nzp, dzt, dn0, a_press, a_scr1, a_scr2,                       &
                     a_rc, a_wp, a_rp, a_rt, a_rh, a_tp,                                 &
                     a_naerop, a_naerot, a_ncloudp, a_ncloudt, a_nprecpp, a_nprecpt,    &
                     a_maerop, a_maerot, a_mcloudp, a_mcloudt, a_mprecpp, a_mprecpt,    &
                     a_nicep,  a_nicet,  a_micep,  a_micet,                             &
                     a_nsnowp, a_nsnowt, a_msnowp, a_msnowt,                            &
                     a_gaerop, a_gaerot, a_rhop,  a_dn,  a_nactd,  a_vactd,   prtcl,    &
                     a_Radry,  a_Rawet,  a_Rcdry,   a_Rcwet,   a_Rpdry,   a_Rpwet,      &
                     a_Ridry,  a_Riwet,  a_Rsdry,   a_Rswet,                            &
                     a_rt,a_rp,a_tt, sst, &
                     a_rsi, a_temp0


    use stat, only : sflg, statistics, acc_massbudged
    use sgsm, only : diffuse
    use srfc, only : surface
    use thrm, only : thermo
    use mcrp, only : micro
    use prss, only : poisson
    use advf, only : fadvect, newdroplet
    use advl, only : ladvect
    use forc, only : forcings
    USE util, ONLY : maskactiv !Juha: Included for SALSA

    USE mo_salsa_driver, ONLY : run_SALSA
    USE mo_submctl, ONLY : nbins,ncld, nprc 
    USE class_ComponentIndex, ONLY : GetNcomp
    USE mpi_interface, ONLY : myid

    logical, intent (out) :: cflflg
    real, intent (out)    :: cflmax

    real :: xtime

    LOGICAL :: zactmask(nzp,nxp,nyp)
    REAL :: zwp(nzp,nxp,nyp), &  !! FOR SINGLE-COLUMN RUNS
            ztkt(nzp,nxp,nyp)
    INTEGER :: zrm
    LOGICAL :: dbg,dbg2,tulostus

    INTEGER :: n4,hhh

    dbg = .FALSE.
    tulostus = .false.
    zwp = 0.5

    xtime = time/86400. + strtim
    cflflg = .false.

    ! The runmode parameter zrm is used by SALSA only
    zrm = 3
    IF ( time < Tspinup ) zrm = 2


    ! Reset ALL tendencies here.
    !----------------------------------------------------------------
    ! "Scalar" timestep
    CALL tend0(.FALSE.)

    ! Put the newly activated to zero
    IF (level >= 4) THEN
       a_vactd = 0.
       a_nactd = 0.
    END IF

    IF (level >= 4 .AND. time < 1.) THEN
       CALL thermo(level)
       CALL SALSA_diagnostics
    END IF

    call surface(sst)

    call diffuse

    call sponge(0)

    if (level >= 1) then

       call thermo(level)

       call forcings(xtime,cntlat,sst)

       IF (level >= 4) THEN

          n4 = GetNcomp(prtcl) + 1 ! Aerosol components + water

          ! Rate of change in absolute temperature (for some ice processes)
          if (time >= 1.) then
             ztkt = a_scr1-a_temp0
             a_temp0 = a_scr1
          else if (time == 0.) then
             a_temp0 = a_scr1
             ztkt = 0.
          end if

          IF ( nxp ==5 .and. nyp == 5 ) THEN
             ! 1D -runs
             CALL run_SALSA(nxp,nyp,nzp,n4,a_press,a_scr1,ztkt,a_rp,a_rt,a_scr2,a_rsi,zwp,a_dn,  &
                  a_naerop,  a_naerot,  a_maerop,  a_maerot,   &
                  a_ncloudp, a_ncloudt, a_mcloudp, a_mcloudt,  &
                  a_nprecpp, a_nprecpt, a_mprecpp, a_mprecpt,  &
                  a_nicep,   a_nicet,   a_micep,   a_micet,    &
                  a_nsnowp,  a_nsnowt,  a_msnowp,  a_msnowt,   &
                  a_nactd,   a_vactd,   a_gaerop,  a_gaerot,   &
                  a_Radry,   a_Rcdry,   a_Rpdry,               &
                  a_Ridry,   a_Rsdry,                          &
                  a_Rawet,   a_Rcwet,   a_Rpwet,               &
                  a_Riwet,   a_Rswet,                          &
                  a_rhop, zrm, prtcl, dtlt, dbg2, time, level  )
          ELSE
             !! for 2D or 3D runs
             CALL run_SALSA(nxp,nyp,nzp,n4,a_press,a_scr1,ztkt,a_rp,a_rt,a_scr2,a_rsi,a_wp,a_dn,  &
                  a_naerop,  a_naerot,  a_maerop,  a_maerot,   &
                  a_ncloudp, a_ncloudt, a_mcloudp, a_mcloudt,  &
                  a_nprecpp, a_nprecpt, a_mprecpp, a_mprecpt,  &
                  a_nicep,   a_nicet,   a_micep,   a_micet,    &
                  a_nsnowp,  a_nsnowt,  a_msnowp,  a_msnowt,   &
                  a_nactd,   a_vactd,   a_gaerop,  a_gaerot,   &
                  a_Radry,   a_Rcdry,   a_Rpdry,               &
                  a_Ridry,   a_Rsdry,                          &
                  a_Rawet,   a_Rcwet,   a_Rpwet,               &
                  a_Riwet,   a_Rswet,                          &
                  a_rhop, zrm, prtcl, dtlt, dbg2, time, level  )
             
          END IF !nxp==5 and nyp == 5
          
       END IF

    end if ! level

    IF (level >= 4)  &
         CALL tend_constrain(n4) ! Tästä pitäisi päästä eroon!

    call update_sclrs

    !-------------------------------------------
    ! "Deposition" timestep
    ! -- Reset only scalar tendencies
    CALL tend0(.TRUE.)

    ! Dont perform sedimentation during spinup for level 4 OR level 5!
    IF (zrm == 3 .OR. level < 4) &
         CALL micro(level)

    IF (level >= 4) CALL tend_constrain(n4)
    CALL update_sclrs

    !-------------------------------------------
    ! "Advection" timestep
    ! -- Reset only scalar tendencies
    call tend0(.TRUE.)

    ! Mask for cloud base activation
    IF (level >= 4)  CALL maskactiv(zactmask,nxp,nyp,nzp,nbins,2,prtcl,a_rh,              &
                                    rc = a_rc,pa_naerop = a_naerop, pa_maerop = a_maerop, &
                                    pt = a_scr1, Rpwet=a_Rawet, w=a_wp, &
                                    pa_ncloud= a_ncloudp(:,:,:,:) )
    ! Get tendencies from cloud base activation
    IF (level >= 4) CALL newdroplet(zactmask)

    CALL fadvect

    IF (level >= 4)  &
         CALL tend_constrain(n4) ! Tästä pitäisi päästä eroon!

    CALL update_sclrs

    CALL thermo(level)

    IF (level >= 4)  &
         CALL SALSA_diagnostics

    call thermo(level)

    ! TÄSTÄ TULEE VÄLILLÄ JOTAIN HÄRÖJÄ + EI OLE TÄRKEÄ? Flagi namelistissa FALSE
    call corlos

    call ladvect

    call buoyancy

    call sponge(1)

    call poisson

    call cfl (cflflg, cflmax)

    CALL thermo(level)

    if (sflg) then
       call statistics (time+dtl)
    end if

    IF (level >= 4) &
       CALL SALSA_diagnostics

  end subroutine t_step
  !
  !----------------------------------------------------------------------
  ! subroutine tend0: sets all tendency arrays to zero
  !
  subroutine tend0(sclonly)

    use grid, only : a_ut, a_vt, a_wt, nscl, a_st, nxyzp, newsclr
    use util, only : azero

    LOGICAL, INTENT(in) :: sclonly ! If true, only put scalar tendencies to zero

    integer :: n

    IF( .NOT. sclonly) &
         call azero(nxyzp,a_ut,a2=a_vt,a3=a_wt)
    do n=1,nscl
       call newsclr(n)
       call azero(nxyzp,a_st)
    end do

  end subroutine tend0
  !
  !----------------------------------------------------------------------
  ! In case of negative tendencies to SALSA arrays, put some constrains
  ! in order to avoid concentrations going negative. This will possibly
  ! slightly affect the conservation of mass - needs testing/revision
  ! Juha Tonttila, FMI, 2014
  !
  SUBROUTINE tend_constrain(nn)

    USE grid, ONLY : a_naerop, a_naerot, a_ncloudp, a_ncloudt, a_nprecpp, a_nprecpt,   &
                     a_maerop, a_maerot, a_mcloudp, a_mcloudt, a_mprecpp, a_mprecpt,   &
                     a_nicep,  a_nicet, a_nsnowp, a_nsnowt,                            & ! ice'n'snow
                     a_micep,  a_micet, a_msnowp, a_msnowt,                            & ! 'ice'n'snow
                     a_gaerop, a_gaerot, dtlt, prtcl,                                   &
                     nxp,nyp,nzp
    USE mo_submctl, ONLY : nbins, ncld, nprc, &
                               nice, nsnw          !ice'n'snow

    INTEGER, INTENT(in) :: nn

    INTEGER :: cc, ii,jj,kk,ni
    REAL :: volp(nn), volt(nn)

    DO jj = 3,nyp-2

       DO ii = 3,nxp-2

          DO kk = 1,nzp

             ! Aerosols
             DO cc = 1,nbins

                !volp = 0.
                !volt = 0.
                !DO ni = 1,nn
                !   volp(ni) = a_maerop(kk,ii,jj,(ni-1)*nbins+cc)
                !   volt(ni) = a_maerop(kk,ii,jj,(ni-1)*nbins+cc)
                !END DO

                IF ( a_naerop(kk,ii,jj,cc)+a_naerot(kk,ii,jj,cc)*dtlt < 0. ) THEN

                   a_naerot(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_naerop(kk,ii,jj,cc))/dtlt,a_naerot(kk,ii,jj,cc))
                   DO ni = 1,nn
                      a_maerot(kk,ii,jj,(ni-1)*nbins+cc) = MAX( ((1.e-10-1.0)*a_maerop(kk,ii,jj,(ni-1)*nbins+cc))/dtlt,  &
                                                               a_maerot(kk,ii,jj,(ni-1)*nbins+cc) )
                   END DO

                END IF

             END DO

             ! Cloud droplets
             DO cc = 1,ncld

                !volp = 0.
                !volt = 0.
                !DO ni = 1,nn
                !   volp(ni) = a_maerop(kk,ii,jj,(ni-1)*nbins+cc)
                !   volt(ni) = a_maerop(kk,ii,jj,(ni-1)*nbins+cc)
                !END DO

                IF ( a_ncloudp(kk,ii,jj,cc)+a_ncloudt(kk,ii,jj,cc)*dtlt < 0. ) THEN

                   a_ncloudt(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_ncloudp(kk,ii,jj,cc))/dtlt,a_ncloudt(kk,ii,jj,cc))
                   DO ni = 1,nn
                      a_mcloudt(kk,ii,jj,(ni-1)*ncld+cc) = MAX( ((1.e-10-1.0)*a_mcloudp(kk,ii,jj,(ni-1)*ncld+cc))/dtlt,  &
                                                               a_mcloudt(kk,ii,jj,(ni-1)*ncld+cc) )
                   END DO

                END IF

             END DO

             ! Precipitation
             DO cc = 1,nprc

                !volp = 0.
                !volt = 0.
                !DO ni = 1,nn
                !   volp(ni) = a_maerop(kk,ii,jj,(ni-1)*nbins+cc)
                !   volt(ni) = a_maerop(kk,ii,jj,(ni-1)*nbins+cc)
                !END DO

                IF ( a_nprecpp(kk,ii,jj,cc)+a_nprecpt(kk,ii,jj,cc)*dtlt < 0. ) THEN

                   a_nprecpt(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_nprecpp(kk,ii,jj,cc))/dtlt,a_nprecpt(kk,ii,jj,cc))
                   DO ni = 1,nn
                      a_mprecpt(kk,ii,jj,(ni-1)*nprc+cc) = MAX( ((1.e-10-1.0)*a_mprecpp(kk,ii,jj,(ni-1)*nprc+cc))/dtlt,  &
                                                               a_mprecpt(kk,ii,jj,(ni-1)*nprc+cc) )
                   END DO

                END IF

             END DO

             ! ice particles
             DO cc = 1,nice

                !volp = 0.
                !volt = 0.
                !DO ni = 1,nn
                !   volp(ni) = a_maerop(kk,ii,jj,(ni-1)*nbins+cc)
                !   volt(ni) = a_maerop(kk,ii,jj,(ni-1)*nbins+cc)
                !END DO

                IF ( a_nicep(kk,ii,jj,cc)+a_nicet(kk,ii,jj,cc)*dtlt < 0. ) THEN

                   a_nicet(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_nicep(kk,ii,jj,cc))/dtlt,a_nicet(kk,ii,jj,cc))
                   DO ni = 1,nn
                      a_micet(kk,ii,jj,(ni-1)*ncld+cc) = MAX( ((1.e-10-1.0)*a_micep(kk,ii,jj,(ni-1)*nice+cc))/dtlt,  &
                                                               a_micet(kk,ii,jj,(ni-1)*nice+cc) )
                   END DO

                END IF

             END DO

             ! Snow
             DO cc = 1,nsnw

                !volp = 0.
                !volt = 0.
                !DO ni = 1,nn
                !   volp(ni) = a_maerop(kk,ii,jj,(ni-1)*nbins+cc)
                !   volt(ni) = a_maerop(kk,ii,jj,(ni-1)*nbins+cc)
                !END DO

                IF ( a_nsnowp(kk,ii,jj,cc)+a_nsnowt(kk,ii,jj,cc)*dtlt < 0. ) THEN

                   a_nsnowt(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_nsnowp(kk,ii,jj,cc))/dtlt,a_nsnowt(kk,ii,jj,cc))
                   DO ni = 1,nn
                      a_msnowt(kk,ii,jj,(ni-1)*nprc+cc) = MAX( ((1.e-10-1.0)*a_msnowp(kk,ii,jj,(ni-1)*nsnw+cc))/dtlt,  &
                                                               a_msnowt(kk,ii,jj,(ni-1)*nsnw+cc) )
                   END DO

                END IF

             END DO

          END DO ! kk

       END DO ! ii

    END DO ! jj



    !WHERE ( a_naerop(:,:,:,:)+a_naerot(:,:,:,:)*dtlt < 0. )
    !   a_naerot(:,:,:,:) = MAX(((1.e-10-1.0)*a_naerop(:,:,:,:))/dtlt,a_naerot(:,:,:,:))
    !ENDWHERE

    !WHERE ( a_ncloudp(:,:,:,:)+a_ncloudt(:,:,:,:)*dtlt < 0. )
    !   a_ncloudt(:,:,:,:) = MAX(((1.e-10-1.0)*a_ncloudp(:,:,:,:))/dtlt,a_ncloudt(:,:,:,:))
    !ENDWHERE

    !WHERE ( a_nprecpp(:,:,:,:)+a_nprecpt(:,:,:,:)*dtlt < 0. )
    !   a_nprecpt(:,:,:,:) = MAX(((1.e-10-1.0)*a_nprecpp(:,:,:,:))/dtlt,a_nprecpt(:,:,:,:))
    !ENDWHERE

    !WHERE ( a_maerop(:,:,:,:)+a_maerot(:,:,:,:)*dtlt < 0. )
    !   a_maerot(:,:,:,:) = MAX(((1.e-10-1.0)*a_maerop(:,:,:,:))/dtlt,a_maerot(:,:,:,:))
    !ENDWHERE

    !WHERE ( a_mcloudp(:,:,:,:)+a_mcloudt(:,:,:,:)*dtlt < 0. )
    !   a_mcloudt(:,:,:,:) = MAX(((1.e-10-1.0)*a_mcloudp(:,:,:,:))/dtlt,a_mcloudt(:,:,:,:))
    !ENDWHERE

    !WHERE ( a_mprecpp(:,:,:,:)+a_mprecpt(:,:,:,:)*dtlt < 0. )
    !   a_mprecpt(:,:,:,:) = MAX(((1.e-10-1.0)*a_mprecpp(:,:,:,:))/dtlt,a_mprecpt(:,:,:,:))
    !ENDWHERE

  END SUBROUTINE tend_constrain
  !
  !----------------------------------------------------------------------
  ! Subroutine cfl: Driver for calling CFL computation subroutine
  !
  subroutine cfl(cflflg,cflmax)

    use grid, only : a_up,a_vp,a_wp,nxp,nyp,nzp,dxi,dyi,dzt,dtlt
    use stat, only : fill_scalar

    logical, intent(out) :: cflflg
    real, intent (out)   :: cflmax
    real, parameter :: cflnum=0.95

    cflmax =  cfll(nzp,nxp,nyp,a_up,a_vp,a_wp,dxi,dyi,dzt,dtlt)

    cflflg = (cflmax > cflnum)
    if (cflflg) print *, 'Warning CFL Violation :', cflmax
    call fill_scalar(1,cflmax)

  end subroutine cfl
  !
  !----------------------------------------------------------------------
  ! Subroutine cfll: Checks CFL criteria, brings down the model if the
  ! maximum thershold is exceeded
  !
  real function cfll(n1,n2,n3,u,v,w,dxi,dyi,dzt,dtlt)

    integer, intent (in) :: n1, n2, n3
    real, dimension (n1,n2,n3), intent (in) :: u, v, w
    real, intent (in)    :: dxi,dyi,dzt(n1),dtlt

    integer :: i, j, k
    cfll=0.
    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
             cfll=max(cfll, dtlt*2.* max(abs(u(k,i,j)*dxi),             &
                  abs(v(k,i,j)*dyi), abs(w(k,i,j)*dzt(k))))
          end do
       end do
    end do

  end function cfll
  !
  !----------------------------------------------------------------------
  ! subroutine update_sclrs:  Updates scalars by applying tendency and
  ! boundary conditions
  !
  subroutine update_sclrs

    use grid, only : a_sp, a_st, a_qp, nscl, nxyzp, nxp, nyp, nzp, dzt, &
         dtlt, newsclr, isgstyp
    use sgsm, only : tkeinit
    use util, only : sclrset

    integer :: n

    do n=1,nscl
       call newsclr(n)
       call update(nzp,nxp,nyp,a_sp,a_st,dtlt)
       call sclrset('mixd',nzp,nxp,nyp,a_sp,dzt)
    end do

    if (isgstyp == 2) then
       call tkeinit(nxyzp,a_qp)
    end if

  end subroutine update_sclrs
  !
  ! ----------------------------------------------------------------------
  ! subroutine update:
  !
  subroutine update(n1,n2,n3,a,fa,dt)

    integer, intent(in)   :: n1, n2, n3
    real, intent (in)     :: fa(n1,n2,n3),dt
    real, intent (in out) :: a(n1,n2,n3)
    integer :: i, j, k

    do j=3,n3-2
       do i=3,n2-2
          do k=2,n1-1
             a(k,i,j) = a(k,i,j) + fa(k,i,j)*dt
          end do
       end do
    end do

  end subroutine update
  !
  ! ----------------------------------------------------------------------
  ! subroutine buoyancy:
  !
  subroutine buoyancy

    use grid, only : a_uc, a_vc, a_wc, a_wt, a_rv, a_rc, a_theta, a_scr1, a_scr3, &
         a_rp, nxp, nyp, nzp, dzm, th00, level, pi1
    use stat, only : sflg, comp_tke
    use util, only : ae1mm
    use thrm, only : update_pi1

    real, dimension (nzp) :: awtbar

    IF (level < 4) THEN
       call boyanc(nzp,nxp,nyp,level,a_wt,a_theta,a_rp,th00,a_scr1,a_rv)
    ELSE IF (level >= 4) THEN
       call boyanc(nzp,nxp,nyp,level,a_wt,a_theta,a_rp,th00,a_scr1,a_rc)
    END IF

    call ae1mm(nzp,nxp,nyp,a_wt,awtbar)
    call update_pi1(nzp,awtbar,pi1)

    if (sflg)  call comp_tke(nzp,nxp,nyp,dzm,th00,a_uc,a_vc,a_wc,a_scr1,a_scr3)

  end subroutine buoyancy
  !
  ! ----------------------------------------------------------------------
  ! subroutine boyanc:
  !
  subroutine boyanc(n1,n2,n3,level,wt,th,rt,th00,scr,rx)

    use defs, only: g, ep2

    integer, intent(in) :: n1,n2,n3,level
    real, intent(in)    :: th00,th(n1,n2,n3),  &
                           rt(n1,n2,n3)  ! This is total water mix rat for level < 4
                                         ! and water vapour mix rat for level = 4
    real, intent(in)    :: rx(n1,n2,n3)  ! This should be water vapour mix rat for level < 4
                                         ! and cloud liquid water mix rat for level = 4 (including rain??)
    real, intent(inout) :: wt(n1,n2,n3)
    real, intent(out)   :: scr(n1,n2,n3)

    integer :: k, i, j
    real :: gover2

    gover2  = 0.5*g

    do j=3,n3-2
       do i=3,n2-2
          if (level >= 2 .and. level < 4) then
             do k=1,n1
                scr(k,i,j)=gover2*((th(k,i,j)*(1.+ep2*rx(k,i,j))-th00)       &
                     /th00-(rt(k,i,j)-rx(k,i,j)))
             end do
          else if (level >= 4) then
             do k=1,n1
                scr(k,i,j)=gover2*((th(k,i,j)*(1.+ep2*rt(k,i,j))-th00)       &
                     /th00-(rx(k,i,j)))
             end do
          else
             do k=1,n1
                scr(k,i,j)=gover2*(th(k,i,j)/th00-1.)
             end do
          end if

          do k=2,n1-2
             wt(k,i,j)=wt(k,i,j)+scr(k,i,j)+scr(k+1,i,j)
          end do
       end do
    end do

  end subroutine boyanc
  !
  ! ----------------------------------------------------------------------
  ! subroutine corlos:  This is the coriolis driver, its purpose is to
  ! from the coriolis accelerations for u and v and add them into the
  ! accumulated tendency arrays of ut and vt.
  !
  subroutine corlos

    use defs, only : omega
    use grid, only : a_uc, a_vc, a_ut, a_vt, nxp, nyp, nzp, u0, v0

    logical, save :: initialized = .False.
    real, save    :: fcor

    integer :: i, j, k

    if (corflg) then
       if (.not.initialized) fcor=2.*omega*sin(cntlat*0.01745329)
       do j=3,nyp-2
          do i=3,nxp-2
             do k=2,nzp
                a_ut(k,i,j)=a_ut(k,i,j) - fcor*(v0(k)-0.25*                   &
                     (a_vc(k,i,j)+a_vc(k,i+1,j)+a_vc(k,i,j-1)+a_vc(k,i+1,j-1)))
                a_vt(k,i,j)=a_vt(k,i,j) + fcor*(u0(k)-0.25*                   &
                     (a_uc(k,i,j)+a_uc(k,i-1,j)+a_uc(k,i,j+1)+a_uc(k,i-1,j+1)))
             end do
          end do
       end do
       initialized = .True.
    end if

  end subroutine corlos
!
! ----------------------------------------------------------------------
! subroutine sponge: does the rayleigh friction for the momentum terms,
! and newtonian damping of thermal term the damping is accumulated with the
! other tendencies
!
  subroutine sponge (isponge)

    use grid, only : u0, v0, a_up, a_vp, a_wp, a_tp, a_ut, a_vt, a_wt, a_tt,&
         nfpt, spng_tfct, spng_wfct, nzp, nxp, nyp, th0, th00

    integer, intent (in) :: isponge

    integer :: i, j, k, kk

    if (maxval(spng_tfct) > epsilon(1.) .and. nfpt > 1) then
       do j=3,nyp-2
          do i=3,nxp-2
             do k=nzp-nfpt,nzp-1
                kk = k+1-(nzp-nfpt)
                if (isponge == 0) then
                   a_tt(k,i,j)=a_tt(k,i,j) - spng_tfct(kk)*                   &
                        (a_tp(k,i,j)-th0(k)+th00)
                else
                   a_ut(k,i,j)=a_ut(k,i,j) - spng_tfct(kk)*(a_up(k,i,j)-u0(k))
                   a_vt(k,i,j)=a_vt(k,i,j) - spng_tfct(kk)*(a_vp(k,i,j)-v0(k))
                   a_wt(k,i,j)=a_wt(k,i,j) - spng_wfct(kk)*(a_wp(k,i,j))
                end if
             end do
          end do
       end do
    end if

  end subroutine sponge

  !
  ! ---------------------------------------------------------------------
  ! SALSA_diagnostics: Update properties for the current timestep:
  !                    E.g. if enough water has evaporated from droplets,
  !                    deplete the cloud droplet bins and move CCN material
  !                    back to the aerosol regime.
  !                    In addition, update the diagnostic scalars for total grid-cell
  !                    liquid water contents.
  !
  ! Juha Tonttila, FMI, 2014
  !

  SUBROUTINE SALSA_diagnostics
    USE grid, ONLY : nxp,nyp,nzp,    &
                     a_naerop,a_maerop,a_ncloudp,a_mcloudp,a_nprecpp,a_mprecpp,      &
                     a_gaerop, a_Radry, a_Rcdry, a_Rpdry, a_Rawet, a_Rcwet, a_Rpwet, &
                     a_rc, a_srp,a_snrp, a_nactd, a_vactd, binMixrat, prtcl, a_rhop, a_dn,  &
                     a_rh, a_scr1, a_ri,a_srs,a_snrs,a_rhi,                                      &
                     a_nicep,a_micep,a_nsnowp,a_msnowp,a_Ridry,a_Riwet,a_Rswet,a_Rsdry !! ice'n'snow
    USE mo_submctl, ONLY : nbins,ncld,nprc,ica,fca,icb,fcb,ira,fra,              &
                               in1a,fn1a,in2a,fn2a,in2b,fn2b,                        &
                               nice,nsnw,iia,fia,iib,fib,isa,fsa,                    & !! ice'n'snow
                               rhosu,rhooc,rhono,rhonh,rhoss,rhowa,rhoic,rhosn,      &  !! Jaakko: rhoic added
                               msu,moc,mno,mnh,mss,mwa,avog,pi6,                     &
                               surfw0,surfi0, rg, nlim, prlim, pi
    USE class_ComponentIndex, ONLY : GetIndex, GetNcomp, IsUsed


    IMPLICIT NONE

    INTEGER :: i,j,k,bc,ba,ww,s,sc,sa,str,end,nc,ino,inh,c,nn
    INTEGER :: cidx(8),aidx(8),pidx(8),iidx(8),sidx(8)

    REAL :: zvol
    REAL :: zvrem
    REAL, PARAMETER :: rth = 1.e-6
    REAL, PARAMETER :: rempty = 1.e-10
    REAL :: zdh2o,zddry
    REAL :: zdiff(fn2a)
    REAL :: ns, bb, aa ! Number of moles, Raoult effect, Kelvin effect; For calculating the critical radius
    REAL :: cdcld(nzp,nxp,nyp,ncld),cdprc(nzp,nxp,nyp,nprc),  & ! Critical diameter for cloud droplets and precipitation
            cdice(nzp,nxp,nyp,nice),cdsnw(nzp,nxp,nyp,nsnw)   ! Critical diameter for cloud droplets and precipitation
    REAL :: cscld(nzp,nxp,nyp,ncld),csprc(nzp,nxp,nyp,nprc)  ! Critical supersaturation for cloud droplets and precip. !!huomhuom
    LOGICAL :: zclosest(fn2a)
    REAL :: vsum

    ! Remove negative values (CHECK HOW MUCH THIS AFFECTS MASS CONCERVATION?)
    a_naerop = MAX(0.,a_naerop)
    a_ncloudp = MAX(0.,a_ncloudp)
    a_nprecpp = MAX(0.,a_nprecpp)
    a_maerop = MAX(0.,a_maerop)
    a_mcloudp = MAX(0.,a_mcloudp)
    a_mprecpp = MAX(0.,a_mprecpp)

    a_nicep = MAX(0.,a_nicep)
    a_nsnowp = MAX(0.,a_nsnowp)
    a_micep = MAX(0.,a_micep)
    a_msnowp = MAX(0.,a_msnowp)

    pidx = 0.
    cidx = 0.
    aidx = 0.
    iidx = 0.
    sidx = 0.
    nn = GetNcomp(prtcl)+1 ! total number of species
    DO s = 1,nn
       pidx(s) = (s-1)*nprc
       cidx(s) = (s-1)*ncld
       aidx(s) = (s-1)*nbins
    END DO

    ! Ghost clouds/precip (DOES THIS AFFECT MASS CONSERVATION?)
    ! Critical radius for cloud and precipitation particles
    DO j = 3,nyp-2
       DO i = 3,nxp-2
          DO k = 1,nzp

             ! Aerosols
             DO c = 1,nbins
                vsum = 0.
                DO s = 1,nn
                   vsum = vsum + a_maerop(k,i,j,(s-1)*nbins+c)
                END DO
                IF (a_naerop(k,i,j,c) > 0. .AND. vsum == 0.) THEN
                   a_naerop(k,i,j,c) = 0.
                   DO s = 1,nn
                      a_maerop(k,i,j,(s-1)*nbins+c) = 0
                   END DO
                END IF
             END DO

             ! Clouds
             DO c = 1,ncld
                vsum = 0.
                DO s = 1,nn
                   vsum = vsum + a_mcloudp(k,i,j,(s-1)*ncld+c)
                END DO
                IF (a_ncloudp(k,i,j,c) > 0. .AND. vsum == 0.) THEN
                   a_ncloudp(k,i,j,c) = 0.
                   DO s = 1,nn
                      a_mcloudp(k,i,j,(s-1)*ncld+c) = 0
                   END DO
                END IF

                ! Critical radius -----------------
                ! Moles of solute
                ns = 0.
                IF (a_ncloudp(k,i,j,c) > prlim) THEN   ! note:use prlim here
                   IF (IsUsed(prtcl,'SO4')) THEN
                      s = GetIndex(prtcl,'SO4')
                      str = (s-1)*ncld + c
                      ns = ns + 3.*a_mcloudp(k,i,j,str)/msu
                   END IF
                   IF (IsUsed(prtcl,'OC')) THEN
                      s = GetIndex(prtcl,'OC')
                      str = (s-1)*ncld + c
                      ns = ns + a_mcloudp(k,i,j,str)/moc
                   END IF
                   IF (IsUsed(prtcl,'NO')) THEN
                      s = GetIndex(prtcl,'NO')
                      str = (s-1)*ncld + c
                      ns = ns + a_mcloudp(k,i,j,str)/mno
                   END IF
                   IF (IsUsed(prtcl,'NH')) THEN
                      s = GetIndex(prtcl,'NH')
                      str = (s-1)*ncld + c
                      ns = ns + a_mcloudp(k,i,j,str)/mnh
                   END IF
                   IF (IsUsed(prtcl,'SS')) THEN
                      s = GetIndex(prtcl,'SS')
                      str = (s-1)*ncld + c
                      ns = ns + 2.*a_mcloudp(k,i,j,str)/mss
                   END IF
                   ns = ns/a_ncloudp(k,i,j,c)

                   bb = 3.*mwa*ns/(4.*pi*rhowa)
                   aa = 4.*mwa*surfw0/(rg*rhowa*a_scr1(k,i,j))
                   cdcld(k,i,j,c) = SQRT(3.*bb/aa)
                ELSE
                   cdcld(k,i,j,c) = rempty
                END IF ! nlim
                ! -----------------------------------

             END DO ! ncld

             ! Precipitation
             DO c = 1,nprc
                IF (a_nprecpp(k,i,j,c) > 0. .AND. a_mprecpp(k,i,j,(nn-1)*nprc+c) == 0.) THEN
                   a_nprecpp(k,i,j,c) = 0.
                   DO s = 1,nn
                      a_mprecpp(k,i,j,(s-1)*nprc+c) = 0.
                   END DO
                END IF

                ! Critical radius -----------------
                ! Moles of solute
                ns = 0.
                IF (a_nprecpp(k,i,j,c) > prlim) THEN
                   IF (IsUsed(prtcl,'SO4')) THEN
                      s = GetIndex(prtcl,'SO4')
                      str = (s-1)*nprc + c
                      ns = ns + 3.*a_mprecpp(k,i,j,str)/msu
                   END IF
                   IF (IsUsed(prtcl,'OC')) THEN
                      s = GetIndex(prtcl,'OC')
                      str = (s-1)*nprc + c
                      ns = ns + a_mprecpp(k,i,j,str)/moc
                   END IF
                   IF (IsUsed(prtcl,'NO')) THEN
                      s = GetIndex(prtcl,'NO')
                      str = (s-1)*nprc + c
                      ns = ns + a_mprecpp(k,i,j,str)/mno
                   END IF
                   IF (IsUsed(prtcl,'NH')) THEN
                      s = GetIndex(prtcl,'NH')
                      str = (s-1)*nprc + c
                      ns = ns + a_mprecpp(k,i,j,str)/mnh
                   END IF
                   IF (IsUsed(prtcl,'SS')) THEN
                      s = GetIndex(prtcl,'SS')
                      str = (s-1)*nprc + c
                      ns = ns + a_mprecpp(k,i,j,str)/mss
                   END IF
                   ns = ns/a_nprecpp(k,i,j,c)

                   bb = 3.*mwa*ns/(4.*pi*rhowa)
                   aa = 4.*mwa*surfw0/(rg*rhowa*a_scr1(k,i,j))
                   cdprc(k,i,j,c) = SQRT(3.*bb/aa)
                ELSE
                   cdprc(k,i,j,c) = rempty
                END IF !prlim
                ! -----------------------------------

             END DO ! nprc

             ! Ice
             DO c = 1,nice
                vsum = 0.
                DO s = 1,nn
                   vsum = vsum + a_micep(k,i,j,(s-1)*nice+c)
                END DO
                IF (a_nicep(k,i,j,c) > 0. .AND. vsum == 0.) THEN
                   a_nicep(k,i,j,c) = 0.
                   DO s = 1,nn
                      a_micep(k,i,j,(s-1)*nice+c) = 0
                   END DO
                END IF

                ! Critical radius -----------------
                ! Moles of solute
                ns = 0.
                IF (a_nicep(k,i,j,c) > nlim) THEN   ! note:use prlim here
                   IF (IsUsed(prtcl,'SO4')) THEN
                      s = GetIndex(prtcl,'SO4')
                      str = (s-1)*nice + c
                      ns = ns + 3.*a_micep(k,i,j,str)/msu
                   END IF
                   IF (IsUsed(prtcl,'OC')) THEN
                      s = GetIndex(prtcl,'OC')
                      str = (s-1)*nice + c
                      ns = ns + a_micep(k,i,j,str)/moc
                   END IF
                   IF (IsUsed(prtcl,'NO')) THEN
                      s = GetIndex(prtcl,'NO')
                      str = (s-1)*nice + c
                      ns = ns + a_micep(k,i,j,str)/mno
                   END IF
                   IF (IsUsed(prtcl,'NH')) THEN
                      s = GetIndex(prtcl,'NH')
                      str = (s-1)*nice + c
                      ns = ns + a_micep(k,i,j,str)/mnh
                   END IF
                   IF (IsUsed(prtcl,'SS')) THEN
                      s = GetIndex(prtcl,'SS')
                      str = (s-1)*nice + c
                      ns = ns + 2.*a_micep(k,i,j,str)/mss
                   END IF
                   ns = ns/a_nicep(k,i,j,c)

                   bb = 3.*mwa*ns/(4.*pi*rhoic)
                   aa = 4.*mwa*surfi0/(rg*rhoic*a_scr1(k,i,j))
                   cdice(k,i,j,c) = max(rempty,SQRT(3.*bb/aa))

                ELSE
                   cdice(k,i,j,c) = rempty
                END IF ! nlim
                ! -----------------------------------

             END DO ! ncld

             ! Snowing
             DO c = 1,nsnw
                IF (a_nsnowp(k,i,j,c) > 0. .AND. a_msnowp(k,i,j,(nn-1)*nsnw+c) == 0.) THEN
                   a_nsnowp(k,i,j,c) = 0.
                   DO s = 1,nn
                      a_msnowp(k,i,j,(s-1)*nsnw+c) = 0.
                   END DO
                END IF

                ! Critical radius -----------------
                ! Moles of solute
                ns = 0.
                IF (a_nsnowp(k,i,j,c) > nlim) THEN
                   IF (IsUsed(prtcl,'SO4')) THEN
                      s = GetIndex(prtcl,'SO4')
                      str = (s-1)*nsnw + c
                      ns = ns + 3.*a_msnowp(k,i,j,str)/msu
                   END IF
                   IF (IsUsed(prtcl,'OC')) THEN
                      s = GetIndex(prtcl,'OC')
                      str = (s-1)*nsnw + c
                      ns = ns + a_msnowp(k,i,j,str)/moc
                   END IF
                   IF (IsUsed(prtcl,'NO')) THEN
                      s = GetIndex(prtcl,'NO')
                      str = (s-1)*nsnw + c
                      ns = ns + a_msnowp(k,i,j,str)/mno
                   END IF
                   IF (IsUsed(prtcl,'NH')) THEN
                      s = GetIndex(prtcl,'NH')
                      str = (s-1)*nsnw + c
                      ns = ns + a_msnowp(k,i,j,str)/mnh
                   END IF
                   IF (IsUsed(prtcl,'SS')) THEN
                      s = GetIndex(prtcl,'SS')
                      str = (s-1)*nsnw + c
                      ns = ns + a_msnowp(k,i,j,str)/mss
                   END IF
                   ns = ns/a_nsnowp(k,i,j,c)

                   bb = 3.*mwa*ns/(4.*pi*rhoic)
                   aa = 4.*mwa*surfi0/(rg*rhoic*a_scr1(k,i,j))
                   cdsnw(k,i,j,c) = max(rempty,SQRT(3.*bb/aa))
                ELSE
                   cdsnw(k,i,j,c) = rempty
                END IF !prlim
                ! -----------------------------------

             END DO ! nsnw

          END DO !k
       END DO !i
    END DO !j

    !
    DO j = 3,nyp-2
       DO i = 3,nxp-2
          DO k = 1,nzp

             !!!!!!!!!!!!!!!!!!!!!!!
             ! CLOUD DROPLETS
             !!!!!!!!!!!!!!!!!!!!!!!

             ! Loop over cloud droplet bins, regime a
             DO bc = ica%cur,fca%cur
                ba = ica%par + (bc-ica%cur) ! Index for parallel aerosol bin
                nc = GetIndex(prtcl,'H2O')
                ww = (nc-1)*ncld + bc ! Index for water content in the cloud bin

                ! Hydrometeor dry radius
                ! If used, remove NO and NH from the total volume because of their volatility
                zvrem = 0
                IF ( IsUsed(prtcl,'NO') ) THEN
                   ino = GetIndex(prtcl,'NO')
                   str = (ino-1)*ncld + bc
                   zvrem = zvrem + a_mcloudp(k,i,j,str)
                END IF
                IF ( IsUsed(prtcl,'NH') ) THEN
                   inh = GetIndex(prtcl,'NH')
                   str = (inh-1)*ncld + bc
                   zvrem = zvrem + a_mcloudp(k,i,j,str)
                END IF

                IF ( a_ncloudp(k,i,j,bc)*a_dn(k,i,j) > prlim ) THEN ! note:prlim

                   CALL binMixrat('cloud','wet',bc,i,j,k,zvol)
                   zvol = zvol/rhowa
                   zdh2o = (zvol/a_ncloudp(k,i,j,bc)/pi6)**(1./3.)

                   ! Loose the droplets if smaller than the critical size
                   IF ( zdh2o < MAX(0.2*cdcld(k,i,j,bc),2.e-6) .AND. a_rh(k,i,j)<0.999 ) THEN

                      ! Move the number of particles from cloud to aerosol bins
                      a_naerop(k,i,j,ba) = a_naerop(k,i,j,ba) + a_ncloudp(k,i,j,bc)
                      a_ncloudp(k,i,j,bc) = 0.

                      ! Move ccn material back to aerosol regime (including water)
                      DO s = 1,GetNcomp(prtcl)+1
                         sc = cidx(s) + bc
                         sa = aidx(s) + ba
                         a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_mcloudp(k,i,j,sc)
                         a_mcloudp(k,i,j,sc) = 0.
                      END DO

                   END IF ! critical radius

                END IF  ! prlim

             END DO ! bc

             ! Loop over cloud droplet bins, regime b
             DO bc = icb%cur,fcb%cur
                ba = icb%par + (bc-icb%cur) ! Index for parallel aerosol bin
                nc = GetIndex(prtcl,'H2O')
                ww = (nc-1)*ncld + bc ! Index for water content in each bin

                ! Hydrometeor dry radius
                ! If used, remove NO and NH from the total volume because of their volatility
                zvrem = 0
                IF ( IsUsed(prtcl,'NO') ) THEN
                   ino = GetIndex(prtcl,'NO')
                   str = (ino-1)*ncld + bc
                   zvrem = zvrem + a_mcloudp(k,i,j,str)
                END IF
                IF ( IsUsed(prtcl,'NH') ) THEN
                   inh = GetIndex(prtcl,'NH')
                   str = (inh-1)*ncld + bc
                   zvrem = zvrem + a_mcloudp(k,i,j,str)
                END IF

                IF ( a_ncloudp(k,i,j,bc)*a_dn(k,i,j) > prlim ) THEN

                   CALL binMixrat('cloud','wet',bc,i,j,k,zvol)
                   zvol = zvol/rhowa
                   zdh2o = (zvol/a_ncloudp(k,i,j,bc)/pi6)**(1./3.)

                   ! Loose the droplets if smaller than critical radius
                   IF ( zdh2o < MAX(0.2*cdcld(k,i,j,bc),2.e-6) .AND. a_rh(k,i,j) < 0.999 ) THEN

                      ! Move the number of particles from cloud to aerosol bins
                      a_naerop(k,i,j,ba) = a_naerop(k,i,j,ba) + a_ncloudp(k,i,j,bc)
                      a_ncloudp(k,i,j,bc) = 0.

                      ! Move ccn material back to aerosol regime (including water)
                      DO s = 1,GetNcomp(prtcl) + 1
                         sc = cidx(s) + bc
                         sa = aidx(s) + ba
                         a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_mcloudp(k,i,j,sc)
                         a_mcloudp(k,i,j,sc) = 0.
                      END DO

                   END IF ! critical radius

                END IF ! prlim

             END DO ! bc

             !!!!!!!!!!!!!!!!!!!!!!!
             ! RAIN DROPS
             !!!!!!!!!!!!!!!!!!!!!!!

             DO bc = ira,fra
                nc = GetIndex(prtcl,'H2O')
                ww = (nc-1)*nprc + bc ! Index for water content in each bin

                ! Hydrometeor dry radius
                ! If used, remove NO and NH from the total volume because of their volatility
                zvrem = 0
                IF ( IsUsed(prtcl,'NO') ) THEN
                   ino = GetIndex(prtcl,'NO')
                   str = (ino-1)*nprc + bc
                   zvrem = zvrem + a_mprecpp(k,i,j,str)
                END IF
                IF ( IsUsed(prtcl,'NH') ) THEN
                   inh = GetIndex(prtcl,'NH')
                   str = (inh-1)*nprc + bc
                   zvrem = zvrem + a_mprecpp(k,i,j,str)
                END IF

                IF ( a_nprecpp(k,i,j,bc)*a_dn(k,i,j) > prlim ) THEN

                   CALL binMixrat('precp','wet',bc,i,j,k,zvol)
                   zvol = zvol/rhowa
                   zdh2o = (zvol/a_nprecpp(k,i,j,bc)/pi6)**(1./3.)

                   ! Loose the droplets if smaller than critical radius
                   IF ( zdh2o < MAX(0.02*cdprc(k,i,j,bc),2.e-6) .AND. a_rh(k,i,j)<0.999 ) THEN

                      ! Move evaporating rain drops to a soluble aerosoln bin with
                      ! the closest match in dry particle radius. Ain't perfect but
                      ! the bin update subroutine in SALSA will take care of the rest.
	                  zvol = 0.
                      zclosest = .FALSE.
                      CALL binMixrat('precp','dry',bc,i,j,k,zvol)
                      zvol = zvol/rhosu
                      zddry = (zvol/a_nprecpp(k,i,j,bc)/pi6)**(1./3.)
                      zdiff(in1a:fn2a) = ABS(2.*a_Radry(k,i,j,in1a:fn2a) - zddry)
                      zclosest(in1a:fn2a) = ( zdiff(in1a:fn2a) == MINVAL(zdiff(in1a:fn2a)) )
                      IF (ALL(zclosest .EQV. .FALSE.)) STOP 'FAIL: zclosest based on NANs'
                      ba = 1
                      DO WHILE( .NOT. zclosest(ba))
                         ba = ba + 1
                      END DO

                      ! Move the number of particles from cloud to aerosol bins
                      a_naerop(k,i,j,ba) = a_naerop(k,i,j,ba) + a_nprecpp(k,i,j,bc)
                      a_nprecpp(k,i,j,bc) = 0.

                      ! Move ccn material back to aerosol regime (including water)
                      DO s = 1,GetNcomp(prtcl) + 1
                         sc = pidx(s) + bc
                         sa = aidx(s) + ba
                         a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_mprecpp(k,i,j,sc)
                         a_mprecpp(k,i,j,sc) = 0.
                      END DO

                   END IF ! Critical radius

                END IF ! prlim

             END DO ! bc

             !!!!!!!!!!!!!!!!!!!!!!!
             ! ICE PARTICLES
             !!!!!!!!!!!!!!!!!!!!!!!

             ! Loop over cloud droplet bins, regime a
             DO bc = iia%cur,fia%cur
                ba = iia%par + (bc-iia%cur) ! Index for parallel aerosol bin
                nc = GetIndex(prtcl,'H2O')
                ww = (nc-1)*nice + bc ! Index for water content in the cloud bin

                ! Hydrometeor dry radius
                ! If used, remove NO and NH from the total volume because of their volatility
                zvrem = 0
                IF ( IsUsed(prtcl,'NO') ) THEN
                   ino = GetIndex(prtcl,'NO')
                   str = (ino-1)*nice + bc
                   zvrem = zvrem + a_micep(k,i,j,str)
                END IF
                IF ( IsUsed(prtcl,'NH') ) THEN
                   inh = GetIndex(prtcl,'NH')
                   str = (inh-1)*nice + bc
                   zvrem = zvrem + a_micep(k,i,j,str)
                END IF

                IF ( a_nicep(k,i,j,bc)*a_dn(k,i,j) > prlim ) THEN ! note:prlim

                   CALL binMixrat('ice','wet',bc,i,j,k,zvol)
				   zvol = zvol/rhoic
                   zdh2o = (zvol/a_nicep(k,i,j,bc)/pi6)**(1./3.)

                   ! Loose the droplets if smaller than the critical size !!huomhuom ice'n'snow
                   IF ( zdh2o < MAX(0.2*cdice(k,i,j,bc),2.e-6) .AND. a_rhi(k,i,j)<0.999 ) THEN

                      ! Move the number of particles from cloud to aerosol bins
                      a_naerop(k,i,j,ba) = a_naerop(k,i,j,ba) + a_nicep(k,i,j,bc)
                      a_nicep(k,i,j,bc) = 0.

                      ! Move ccn material back to aerosol regime (including water)
                      DO s = 1,GetNcomp(prtcl)+1
                         sc = iidx(s) + bc
                         sa = aidx(s) + ba
                         a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_micep(k,i,j,sc)
                         a_micep(k,i,j,sc) = 0.
                      END DO

                   END IF ! critical radius

                END IF  ! prlim

             END DO ! bc

             ! Loop over cloud droplet bins, regime b
             DO bc = iib%cur,fib%cur
                ba = iib%par + (bc-iib%cur) ! Index for parallel aerosol bin
                nc = GetIndex(prtcl,'H2O')
                ww = (nc-1)*nice + bc ! Index for water content in each bin

                ! Hydrometeor dry radius
                ! If used, remove NO and NH from the total volume because of their volatility
                zvrem = 0
                IF ( IsUsed(prtcl,'NO') ) THEN
                   ino = GetIndex(prtcl,'NO')
                   str = (ino-1)*nice + bc
                   zvrem = zvrem + a_micep(k,i,j,str)
                END IF
                IF ( IsUsed(prtcl,'NH') ) THEN
                   inh = GetIndex(prtcl,'NH')
                   str = (inh-1)*nice + bc
                   zvrem = zvrem + a_micep(k,i,j,str)
                END IF

                IF ( a_nicep(k,i,j,bc)*a_dn(k,i,j) > nlim ) THEN

                   CALL binMixrat('ice','wet',bc,i,j,k,zvol)
				   zvol = zvol/rhoic
                   zdh2o = (zvol/a_nicep(k,i,j,bc)/pi6)**(1./3.)

                   ! Loose the droplets if smaller than critical radius !!huomhuom ice'n'snow
                   IF ( zdh2o < MAX(0.2*cdice(k,i,j,bc),2.e-6) .AND. a_rhi(k,i,j) < 0.999 ) THEN

                      ! Move the number of particles from cloud to aerosol bins
                      a_naerop(k,i,j,ba) = a_naerop(k,i,j,ba) + a_nicep(k,i,j,bc)
                      a_nicep(k,i,j,bc) = 0.

                      ! Move ccn material back to aerosol regime (including water)
                      DO s = 1,GetNcomp(prtcl) + 1
                         sc = iidx(s) + bc
                         sa = aidx(s) + ba
                         a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_micep(k,i,j,sc)
                         a_micep(k,i,j,sc) = 0.
                      END DO

                   END IF ! critical radius

                END IF ! nlim

             END DO ! bc

             !!!!!!!!!!!!!!!!!!!!!!!
             ! SNOW FLAKES
             !!!!!!!!!!!!!!!!!!!!!!!

             DO bc = isa,fsa
                nc = GetIndex(prtcl,'H2O')
                ww = (nc-1)*nsnw + bc ! Index for water content in each bin

                ! Hydrometeor dry radius
                ! If used, remove NO and NH from the total volume because of their volatility
                zvrem = 0
                IF ( IsUsed(prtcl,'NO') ) THEN
                   ino = GetIndex(prtcl,'NO')
                   str = (ino-1)*nsnw + bc
                   zvrem = zvrem + a_msnowp(k,i,j,str)
                END IF
                IF ( IsUsed(prtcl,'NH') ) THEN
                   inh = GetIndex(prtcl,'NH')
                   str = (inh-1)*nsnw + bc
                   zvrem = zvrem + a_msnowp(k,i,j,str)
                END IF

                IF ( a_nsnowp(k,i,j,bc)*a_dn(k,i,j) > nlim ) THEN

                   CALL binMixrat('snow','wet',bc,i,j,k,zvol)
				   zvol = zvol/rhosn
                   zdh2o = (zvol/a_nsnowp(k,i,j,bc)/pi6)**(1./3.)

                   ! Loose the droplets if smaller than critical radius !!huomhuom a_rhi ice'n'snow
                   IF ( zdh2o < MAX(0.02*cdsnw(k,i,j,bc),2.e-6) .AND. a_rhi(k,i,j)<0.999 ) THEN

                      ! Move evaporating rain drops to a soluble aerosoln bin with
                      ! the closest match in dry particle radius. Ain't perfect but
                      ! the bin update subroutine in SALSA will take care of the rest.
                      zclosest = .FALSE.
                      CALL binMixrat('snow','dry',bc,i,j,k,zvol)
					  zvol = zvol/rhosn
                      zddry = (zvol/a_nsnowp(k,i,j,bc)/pi6)**(1./3.)
                      zdiff(in1a:fn2a) = ABS(2.*a_Radry(k,i,j,in1a:fn2a) - zddry)
                      zclosest(in1a:fn2a) = ( zdiff(in1a:fn2a) == MINVAL(zdiff(in1a:fn2a)) )
                      IF (ALL(zclosest .EQV. .FALSE.)) STOP 'FAIL: zclosest based on NANs'
                      ba = 1
                      DO WHILE( .NOT. zclosest(ba))
                         ba = ba + 1
                      END DO

                      ! Move the number of particles from cloud to aerosol bins
                      a_naerop(k,i,j,ba) = a_naerop(k,i,j,ba) + a_nsnowp(k,i,j,bc)
                      a_nsnowp(k,i,j,bc) = 0.

                      ! Move ccn material back to aerosol regime (including water)
                      DO s = 1,GetNcomp(prtcl) + 1
                         sc = sidx(s) + bc
                         sa = aidx(s) + ba
                         a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_msnowp(k,i,j,sc)
                         a_msnowp(k,i,j,sc) = 0.
                      END DO

                   END IF ! Critical radius

                END IF ! nlim

             END DO ! bc

             !!!!!!!!!!!!!!!!!!!!!!
             ! AEROSOLS
             !!!!!!!!!!!!!!!!!!!!!!
             ! MIHIS VESI HIUKKASTEN POISTUESSA??

             ! Loop over aerosol bins, regime 1a
             DO ba = in1a,fn1a

                zvrem = 0
                IF ( IsUsed(prtcl,'NO') ) THEN
                   ino = GetIndex(prtcl,'NO')
                   str = (ino-1)*nbins + ba
                   zvrem = zvrem + a_maerop(k,i,j,str)
                END IF
                IF ( IsUsed(prtcl,'NH') ) THEN
                   inh = GetIndex(prtcl,'NH')
                   str = (inh-1)*nbins + ba
                   zvrem = zvrem + a_maerop(k,i,j,str)
                END IF


                IF ( a_naerop(k,i,j,ba) > nlim ) THEN

                   CALL binMixrat('aerosol','dry',ba,i,j,k,zvol)
                   zvol = zvol/rhosu
                   zddry = (zvol/a_naerop(k,i,j,ba)/pi6)**(1./3.)

                   IF ( zddry < 1.e-10  ) THEN

                      IF (IsUsed(prtcl,'SO4')) THEN
                         nc = GetIndex(prtcl,'SO4')
                         s = (nc-1)*nbins + ba
                         a_gaerop(k,i,j,1) = a_gaerop(k,i,j,1) + a_maerop(k,i,j,s) / msu * avog
                         a_maerop(k,i,j,s) = 0.
                      END IF

                      IF (IsUsed(prtcl,'OC')) THEN
                         nc = GetIndex(prtcl,'OC')
                         s = (nc-1)*nbins + ba
                         a_gaerop(k,i,j,5) = a_gaerop(k,i,j,5) + a_maerop(k,i,j,s) / moc * avog
                         a_maerop(k,i,j,s) = 0.
                      END IF

                      IF (IsUsed(prtcl,'NO')) THEN
                         nc = GetIndex(prtcl,'NO')
                         s = (nc-1)*nbins + ba
                         a_gaerop(k,i,j,2) = a_gaerop(k,i,j,2) + a_maerop(k,i,j,s) / mno * avog
                         a_maerop(k,i,j,s) = 0.
                      END IF

                      IF (IsUsed(prtcl,'NH')) THEN
                         nc = GetIndex(prtcl,'NH')
                         s = (nc-1)*nbins + ba
                         a_gaerop(k,i,j,3) = a_gaerop(k,i,j,3) + a_maerop(k,i,j,s) / mnh * avog
                         a_maerop(k,i,j,s) = 0.
                      END IF
                      a_naerop(k,i,j,ba) = 0.

                   END IF ! Rdry

                END IF ! nlim

             END DO ! ba

             ! Loop over aerosol bins, regime 2a
             DO ba = in2a,fn2a

                zvrem = 0
                IF ( IsUsed(prtcl,'NO') ) THEN
                   ino = GetIndex(prtcl,'NO')
                   str = (ino-1)*nbins + ba
                   zvrem = zvrem + a_maerop(k,i,j,str)
                END IF
                IF ( IsUsed(prtcl,'NH') ) THEN
                   inh = GetIndex(prtcl,'NH')
                   str = (inh-1)*nbins + ba
                   zvrem = zvrem + a_maerop(k,i,j,str)
                END IF

                IF ( a_naerop(k,i,j,ba) > nlim ) THEN

                   CALL binMixrat('aerosol','dry',ba,i,j,k,zvol)
                   zvol = zvol/rhosu
                   zddry = (zvol/a_naerop(k,i,j,ba)/pi6)**(1./3.)

                   IF ( zddry < 1.e-10  ) THEN

                      IF (IsUsed(prtcl,'SO4')) THEN
                         nc = GetIndex(prtcl,'SO4')
                         s = (nc-1)*nbins + ba
                         a_gaerop(k,i,j,1) = a_gaerop(k,i,j,1) + a_maerop(k,i,j,s) / msu * avog
                         a_maerop(k,i,j,s) = 0.
                      END IF

                      IF (IsUsed(prtcl,'OC')) THEN
                         nc = GetIndex(prtcl,'OC')
                         s = (nc-1)*nbins + ba
                         a_gaerop(k,i,j,5) = a_gaerop(k,i,j,5) + a_maerop(k,i,j,s) / moc * avog
                         a_maerop(k,i,j,s) = 0.
                      END IF

                      IF (IsUsed(prtcl,'NO')) THEN
                         nc = GetIndex(prtcl,'NO')
                         s = (nc-1)*nbins + ba
                         a_gaerop(k,i,j,2) = a_gaerop(k,i,j,2) + a_maerop(k,i,j,s) / mno * avog
                         a_maerop(k,i,j,s) = 0.
                      END IF

                      IF (IsUsed(prtcl,'NH')) THEN
                         nc = GetIndex(prtcl,'NH')
                         s = (nc-1)*nbins + ba
                         a_gaerop(k,i,j,3) = a_gaerop(k,i,j,3) + a_maerop(k,i,j,s) / mnh * avog
                         a_maerop(k,i,j,s) = 0.
                      END IF
                      a_naerop(k,i,j,ba) = 0.

                   END IF  ! Rdry

                END IF ! nlim

             END DO !ba

             ! Loop over aerosol bins, regime 2b
             DO ba = in2b,fn2b

                zvrem = 0
                IF ( IsUsed(prtcl,'NO') ) THEN
                   ino = GetIndex(prtcl,'NO')
                   str = (ino-1)*nbins + ba
                   zvrem = zvrem + a_maerop(k,i,j,str)
                END IF
                IF ( IsUsed(prtcl,'NH') ) THEN
                   inh = GetIndex(prtcl,'NH')
                   str = (inh-1)*nbins + ba
                   zvrem = zvrem + a_maerop(k,i,j,str)
                END IF

                IF ( a_naerop(k,i,j,ba) > nlim ) THEN

                   CALL binMixrat('aerosol','dry',ba,i,j,k,zvol)
                   zvol = zvol/rhosu
                   zddry = (zvol/a_naerop(k,i,j,ba)/pi6)**(1./3.)

                   IF ( zddry < 1.e-10  ) THEN

                      IF (IsUsed(prtcl,'SO4')) THEN
                         nc = GetIndex(prtcl,'SO4')
                         s = (nc-1)*nbins + ba
                         a_gaerop(k,i,j,1) = a_gaerop(k,i,j,1) + a_maerop(k,i,j,s) / msu * avog
                         a_maerop(k,i,j,s) = 0.
                      END IF

                      IF (IsUsed(prtcl,'OC')) THEN
                         nc = GetIndex(prtcl,'OC')
                         s = (nc-1)*nbins + ba
                         a_gaerop(k,i,j,5) = a_gaerop(k,i,j,5) + a_maerop(k,i,j,s) / moc * avog
                         a_maerop(k,i,j,s) = 0.
                      END IF

                      IF (IsUsed(prtcl,'NO')) THEN
                         nc = GetIndex(prtcl,'NO')
                         s = (nc-1)*nbins + ba
                         a_gaerop(k,i,j,2) = a_gaerop(k,i,j,2) + a_maerop(k,i,j,s) / mno * avog
                         a_maerop(k,i,j,s) = 0.
                      END IF

                      IF (IsUsed(prtcl,'NH')) THEN
                         nc = GetIndex(prtcl,'NH')
                         s = (nc-1)*nbins + ba
                         a_gaerop(k,i,j,3) = a_gaerop(k,i,j,3) + a_maerop(k,i,j,s) / mnh * avog
                         a_maerop(k,i,j,s) = 0.
                      END IF
                      a_naerop(k,i,j,ba) = 0.

                   END IF ! Rdry

                END IF ! nlim

             END DO ! ba

             ! Update particle radiae
             ! Aerosols:
             DO ba = 1,nbins
                IF (a_naerop(k,i,j,ba) > nlim) THEN
                   CALL binMixrat('aerosol','dry',ba,i,j,k,zvol)
                   zvol = zvol/rhosu
                   a_Radry(k,i,j,ba) = 0.5*(zvol/a_naerop(k,i,j,ba)/pi6)**(1./3.)
                   CALL binMixrat('aerosol','wet',ba,i,j,k,zvol)
                   zvol = zvol/1500.
                   a_Rawet(k,i,j,ba) = 0.5*(zvol/a_naerop(k,i,j,ba)/pi6)**(1./3.)
                ELSE
                   a_Radry(k,i,j,ba) = rempty
                   a_Rawet(k,i,j,ba) = rempty
                END IF
             END DO
             DO bc = 1,ncld
                IF (a_ncloudp(k,i,j,bc) > prlim) THEN
                   CALL binMixrat('cloud','dry',bc,i,j,k,zvol)
                   zvol = zvol/rhosu
                   a_Rcdry(k,i,j,bc) = 0.5*(zvol/a_ncloudp(k,i,j,bc)/pi6)**(1./3.)
                   CALL binMixrat('cloud','wet',bc,i,j,k,zvol)
                   zvol = zvol/rhowa
                   a_Rcwet(k,i,j,bc) = 0.5*(zvol/a_ncloudp(k,i,j,bc)/pi6)**(1./3.)
                ELSE
                   a_Rcdry(k,i,j,bc) = rempty
                   a_Rcwet(k,i,j,bc) = rempty
                END IF
             END DO
             DO bc = 1,nprc
                IF (a_nprecpp(k,i,j,bc) > prlim) THEN
                   CALL binMixrat('precp','dry',bc,i,j,k,zvol)
                   zvol = zvol/rhosu
                   a_Rpdry(k,i,j,bc) = 0.5*(zvol/a_nprecpp(k,i,j,bc)/pi6)**(1./3.)
                   CALL binMixrat('precp','wet',bc,i,j,k,zvol)
                   zvol = zvol/rhowa
                   a_Rpwet(k,i,j,bc) = 0.5*(zvol/a_nprecpp(k,i,j,bc)/pi6)**(1./3.)
                   !if (a_Rpwet(k,i,j,bc) > 1  &
                   !½ ) then
                   ! write(*,*) 'step sade märkäsäde aika iso debugkebab', a_Rpwet(k,i,j,bc)
                   ! stop
                   !end if
                ELSE
                   a_Rpdry(k,i,j,bc) = rempty
                   a_Rpwet(k,i,j,bc) = rempty
                END IF
             END DO
             DO bc = 1,nice
                IF (a_nicep(k,i,j,bc) > nlim) THEN
                   CALL binMixrat('ice','dry',bc,i,j,k,zvol)
				   zvol = zvol/rhosu
                   a_Ridry(k,i,j,bc) = 0.5*(zvol/a_nicep(k,i,j,bc)/pi6)**(1./3.)
                   CALL binMixrat('ice','wet',bc,i,j,k,zvol)
				   zvol = zvol/rhoic
                   a_Riwet(k,i,j,bc) = 0.5*(zvol/a_nicep(k,i,j,bc)/pi6)**(1./3.)
                    if (a_Riwet(k,i,j,bc) > 1  &
                    ) then
                    write(*,*) 'step jää märkäsäde aika iso debugkebab'
                    stop
                    end if
                ELSE
                   a_Ridry(k,i,j,bc) = rempty
                   a_Riwet(k,i,j,bc) = rempty
                END IF
             END DO
             DO bc = 1,nsnw
                IF (a_nsnowp(k,i,j,bc) > nlim) THEN
                   CALL binMixrat('snow','dry',bc,i,j,k,zvol)
				   zvol = zvol/rhosu
                   a_Rsdry(k,i,j,bc) = 0.5*(zvol/a_nsnowp(k,i,j,bc)/pi6)**(1./3.)
                   CALL binMixrat('snow','wet',bc,i,j,k,zvol)
				   zvol = zvol/rhosn
                   a_Rswet(k,i,j,bc) = 0.5*(zvol/a_nsnowp(k,i,j,bc)/pi6)**(1./3.)
                ELSE
                   a_Rsdry(k,i,j,bc) = rempty
                   a_Rswet(k,i,j,bc) = rempty
                END IF
             END DO
             ! Update diagnostic tracers
             ! Liquid water content
             nc = GetIndex(prtcl,'H2O')
             ! Aerosols, regime a
             str = (nc-1)*nbins + in1a
             end = (nc-1)*nbins + fn2a
             a_rc(k,i,j) = SUM(a_maerop(k,i,j,str:end))
             ! Aerosols, regime b
             str = (nc-1)*nbins + in2b
             end = (nc-1)*nbins + fn2b
             a_rc(k,i,j) = a_rc(k,i,j) + SUM(a_maerop(k,i,j,str:end))
             ! Clouds, regime a
             str = (nc-1)*ncld+ica%cur
             end = (nc-1)*ncld+fca%cur
             a_rc(k,i,j) = a_rc(k,i,j) + SUM(a_mcloudp(k,i,j,str:end))
             ! Clouds, regime b
             str = (nc-1)*ncld+icb%cur
             end = (nc-1)*ncld+fcb%cur
             a_rc(k,i,j) = a_rc(k,i,j) + SUM(a_mcloudp(k,i,j,str:end))
             ! Precipitating water
             ! regime a
             str = (nc-1)*nprc+ira
             end = (nc-1)*nprc+fra
             a_srp(k,i,j) = SUM(a_mprecpp(k,i,j,str:end))
             a_snrp(k,i,j) = SUM(a_nprecpp(k,i,j,ira:fra))
             a_rhop(k,i,j,:) = 1500.

             ! ice, regime a
             str = (nc-1)*nice+iia%cur
             end = (nc-1)*nice+fia%cur
             a_ri(k,i,j) = a_ri(k,i,j) + SUM(a_micep(k,i,j,str:end))
             ! ice, regime b
             str = (nc-1)*nice+iib%cur
             end = (nc-1)*nice+fib%cur
             a_ri(k,i,j) = a_ri(k,i,j) + SUM(a_micep(k,i,j,str:end))
             ! Precipitating ice  = snow
             ! regime a
             str = (nc-1)*nsnw+isa
             end = (nc-1)*nsnw+fsa
             a_srs(k,i,j) = SUM(a_msnowp(k,i,j,str:end))
	         a_snrs(k,i,j) = SUM(a_nsnowp(k,i,j,isa:fsa))

          END DO   ! k
       END DO   ! i
    END DO   ! j


  END SUBROUTINE SALSA_diagnostics


end module step