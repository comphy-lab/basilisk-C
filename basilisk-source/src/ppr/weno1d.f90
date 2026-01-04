
    !
    ! This program may be freely redistributed under the 
    ! condition that the copyright notices (including this 
    ! entire header) are not removed, and no compensation 
    ! is received through use of the software.  Private, 
    ! research, and institutional use is free.  You may 
    ! distribute modified versions of this code UNDER THE 
    ! CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE 
    ! TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE 
    ! ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE 
    ! MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR 
    ! NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution 
    ! of this code as part of a commercial system is 
    ! permissible ONLY BY DIRECT ARRANGEMENT WITH THE 
    ! AUTHOR.  (If you are not directly supplying this 
    ! code to a customer, and you are instead telling them 
    ! how they can obtain it for free, then you are not 
    ! required to make any arrangement with me.) 
    !
    ! Disclaimer:  Neither I nor: Columbia University, the 
    ! National Aeronautics and Space Administration, nor 
    ! the Massachusetts Institute of Technology warrant 
    ! or certify this code in any way whatsoever.  This 
    ! code is provided "as-is" to be used at your own risk.
    !
    !

    !    
    ! WENO1D.f90: WENO-style slope-limiting for 1d reconst.
    !
    ! Darren Engwirda 
    ! 08-Sep-2016
    ! de2363 [at] columbia [dot] edu
    !
    !

    pure subroutine wenoi (npos,delx,oscl,ipos, &
        &                  ivar,halo,&
        &                  wlim,wval )

    !
    ! NPOS  no. edges over grid.
    ! DELX  grid-cell spacing array. SIZE(DELX) == +1 if
    !       the grid is uniformly spaced .
    ! OSCL  cell-centred oscillation-detectors, where OSCL 
    !       has SIZE = +2-by-NVAR-by-NPOS-1. OSCL is given
    !       by calls to OSCLI().
    ! IPOS  grid-cell index for which to calc. weights .
    ! IVAR  state-var index for which to calc/ weights .
    ! HALO  width of recon. stencil, symmetric about IPOS .
    ! WLIM  limiter treatment at endpoints, monotonic or
    !       otherwise . 
    ! WVAL  WENO weights vector, such that FHAT = WVAL(1) *
    !       UHAT + WVAL(2) * LHAT, where UHAT and LHAT are 
    !       the unlimited and monotonic grid-cell profiles
    !       respectively .
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)  :: npos,halo
        integer, intent(in)  :: ipos,ivar
        integer, intent(in)  :: wlim
        real*8 , intent(in)  :: delx(:)
        real*8 , intent(in)  :: oscl(:,:,:)
        real*8 , intent(out) :: wval(2)
    
    !------------------------------------------- variables !
        real*8 :: omin,omax,wsum
        
        real*8 , parameter :: ZERO = +1.d-16
       
        if (size(delx).gt.+1) then
        
    !------------------- use variable grid spacing variant !
        
        call wenov(npos,delx,oscl, &
        &          ipos,ivar,halo, &
        &          wlim,omin,omax)
        
        else
        
    !------------------- use constant grid spacing variant !
        
        call wenoc(npos,delx,oscl, &
        &          ipos,ivar,halo, &
        &          wlim,omin,omax)
        
        end if

    !------------------ compute WENO-style profile weights !

        omax = omax + ZERO
        omin = omin + ZERO

        if (halo .ge. +3) then

        wval(1) = +1.0d+7 / omax ** 3
        wval(2) = +1.0d+0 / omin ** 3

        else &
    &   if (halo .le. +2) then
    
        wval(1) = +1.0d+5 / omax ** 3
        wval(2) = +1.0d+0 / omin ** 3
    
        end if

        wsum = wval(1) + wval(2) + ZERO
        wval(1) = wval(1) / wsum
    !   wval(2) = wval(2) / wsum
        wval(2) =-wval(1) + 1.d0 ! wval(2)/wsum but robust !

        return

    end  subroutine

    pure subroutine wenov (npos,delx,oscl,ipos, &
        &                  ivar,halo,&
        &                  wlim,omin,omax)

    !
    ! *this is the variable grid-spacing variant .
    !
    ! NPOS  no. edges over grid.
    ! DELX  grid-cell spacing array. SIZE(DELX) == +1 if
    !       the grid is uniformly spaced .
    ! OSCL  cell-centred oscillation-detectors, where OSCL 
    !       has SIZE = +2-by-NVAR-by-NPOS-1. OSCL is given
    !       by calls to OSCLI().
    ! IPOS  grid-cell index for which to calc. weights .
    ! IVAR  state-var index for which to calc/ weights .
    ! HALO  width of recon. stencil, symmetric about IPOS .
    ! WLIM  limiter treatment at endpoints, monotonic or
    !       otherwise . 
    ! OMIN  min. and max. oscillation indicators over the
    ! OMAX  local re-con. stencil .
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)  :: npos,halo
        integer, intent(in)  :: ipos,ivar
        integer, intent(in)  :: wlim
        real*8 , intent(in)  :: delx(:)
        real*8 , intent(in)  :: oscl(:,:,:)
        real*8 , intent(out) :: omin,omax

    !------------------------------------------- variables !
        integer :: hpos
        integer :: head,tail
        integer :: imin,imax
        real*8  :: deli,delh
        real*8  :: hh00,hsqr
        real*8  :: dfx1,dfx2
        real*8  :: oval
        
    !------------------- calc. lower//upper stencil bounds !    

        head = 1; tail = npos - 1

        if(wlim.eq.mono_limit) then

    !---------------------- deactivate WENO at boundaries !
        
        if (ipos-halo.lt.head) then
        
            omax = 1.d0 
            omin = 0.d0 ; return
        
        end if
        
        if (ipos+halo.gt.tail) then
            
            omax = 1.d0 
            omin = 0.d0 ; return
        
        end if

        end if

    !---------------------- truncate stencil at boundaries !
    
       imin = max(ipos-halo,head)
       imax = min(ipos+halo,tail)
      
    !------------------ find min/max indicators on stencil !
 
        dfx1 = oscl(1,ivar,ipos)
        dfx2 = oscl(2,ivar,ipos)
      
        hh00 = delx(ipos+0)**1
        hsqr = delx(ipos+0)**2
      
        oval =(hh00 * dfx1)**2 &
        &    +(hsqr * dfx2)**2
      
        omin = oval
        omax = oval
      
    !---------------------------------------- "lower" part !
      
        delh = 0.d0

        do  hpos = ipos-1, imin, -1
        
    !------------------ calc. derivatives centred on IPOS. !     
   
        deli = delx(hpos+0) &
    &        + delx(hpos+1)
    
        delh = delh + deli*.5d0
        
        dfx1 = oscl(1,ivar,hpos)
        dfx2 = oscl(2,ivar,hpos)

        dfx1 = dfx1 + dfx2*delh

    !------------------ indicator: NORM(H^N * D^N/DX^N(F)) !

        oval = (hh00 * dfx1)**2 &
    &        + (hsqr * dfx2)**2

        if (oval .lt. omin) then
            omin = oval 
        else &
    &   if (oval .gt. omax) then
            omax = oval
        end if
        
        end do
      
    !---------------------------------------- "upper" part !     
      
        delh = 0.d0
        
        do  hpos = ipos+1, imax, +1
        
    !------------------ calc. derivatives centred on IPOS. !     
   
        deli = delx(hpos+0) &
    &        + delx(hpos-1)
    
        delh = delh - deli*.5d0
        
        dfx1 = oscl(1,ivar,hpos)
        dfx2 = oscl(2,ivar,hpos)

        dfx1 = dfx1 + dfx2*delh

    !------------------ indicator: NORM(H^N * D^N/DX^N(F)) !

        oval = (hh00 * dfx1)**2 &
    &        + (hsqr * dfx2)**2

        if (oval .lt. omin) then
            omin = oval 
        else &
    &   if (oval .gt. omax) then
            omax = oval
        end if
        
        end do

        return

    end  subroutine
    
    pure subroutine wenoc (npos,delx,oscl,ipos, &
        &                  ivar,halo,&
        &                  wlim,omin,omax)

    !
    ! *this is the constant grid-spacing variant .
    !
    ! NPOS  no. edges over grid.
    ! DELX  grid-cell spacing array. SIZE(DELX) == +1 if
    !       the grid is uniformly spaced .
    ! OSCL  cell-centred oscillation-detectors, where OSCL 
    !       has SIZE = +2-by-NVAR-by-NPOS-1. OSCL is given
    !       by calls to OSCLI().
    ! IPOS  grid-cell index for which to calc. weights .
    ! IVAR  state-var index for which to calc/ weights .
    ! HALO  width of recon. stencil, symmetric about IPOS .
    ! WLIM  limiter treatment at endpoints, monotonic or
    !       otherwise . 
    ! OMIN  min. and max. oscillation indicators over the
    ! OMAX  local re-con. stencil .
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)  :: npos,halo
        integer, intent(in)  :: ipos,ivar
        integer, intent(in)  :: wlim
        real*8 , intent(in)  :: delx(1)
        real*8 , intent(in)  :: oscl(:,:,:)
        real*8 , intent(out) :: omin,omax

    !------------------------------------------- variables !
        integer :: hpos
        integer :: head,tail
        integer :: imin,imax
        real*8  :: delh
        real*8  :: dfx1,dfx2
        real*8  :: oval
    
    !------------------- calc. lower//upper stencil bounds !    

        head = 1; tail = npos - 1

        if(wlim.eq.mono_limit) then

    !---------------------- deactivate WENO at boundaries !
        
        if (ipos-halo.lt.head) then
        
            omax = 1.d0 
            omin = 0.d0 ; return
        
        end if
        
        if (ipos+halo.gt.tail) then
            
            omax = 1.d0 
            omin = 0.d0 ; return
        
        end if

        end if

    !---------------------- truncate stencil at boundaries !
        
        imin = max(ipos-halo,head)
        imax = min(ipos+halo,tail)
      
    !------------------ find min/max indicators on stencil !
 
        dfx1 = oscl(1,ivar,ipos)
        dfx2 = oscl(2,ivar,ipos)
      
        oval = (2.d0**1*dfx1)**2 &
        &    + (2.d0**2*dfx2)**2
      
        omin = oval
        omax = oval
      
    !---------------------------------------- "lower" part !
      
        delh = 0.d0

        do  hpos = ipos-1, imin, -1
        
    !------------------ calc. derivatives centred on IPOS. !     

        delh = delh + 2.d0
        
        dfx1 = oscl(1,ivar,hpos)
        dfx2 = oscl(2,ivar,hpos)

        dfx1 = dfx1 + dfx2*delh

    !------------------ indicator: NORM(H^N * D^N/DX^N(F)) !

        oval = (2.d0**1*dfx1)**2 &
    &        + (2.d0**2*dfx2)**2

        if (oval .lt. omin) then
            omin = oval 
        else &
    &   if (oval .gt. omax) then
            omax = oval
        end if
        
        end do
      
    !---------------------------------------- "upper" part !     
      
        delh = 0.d0
        
        do  hpos = ipos+1, imax, +1
        
    !------------------ calc. derivatives centred on IPOS. !     

        delh = delh - 2.d0
        
        dfx1 = oscl(1,ivar,hpos)
        dfx2 = oscl(2,ivar,hpos)

        dfx1 = dfx1 + dfx2*delh
    
    !------------------ indicator: NORM(H^N * D^N/DX^N(F)) !

        oval = (2.d0**1*dfx1)**2 &
    &        + (2.d0**2*dfx2)**2
  
        if (oval .lt. omin) then
            omin = oval 
        else &
    &   if (oval .gt. omax) then
            omax = oval
        end if

        end do

        return

    end  subroutine
    
    
    
