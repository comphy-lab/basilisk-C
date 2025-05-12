
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
    ! PPM.f90: 1d slope-limited, piecewise parabolic recon.
    !
    ! Darren Engwirda 
    ! 08-Sep-2016
    ! de2363 [at] columbia [dot] edu
    !
    !
  
    ! P. Colella and PR. Woodward, The Piecewise Parabolic 
    ! Method (PPM) for gas-dynamical simulations., J. Comp. 
    ! Phys., 54 (1), 1984, 174-201,
    ! https://doi.org/10.1016/0021-9991(84)90143-8
    !

    pure subroutine ppm(npos,nvar,ndof,delx, &
        &               fdat,fhat,edge,oscl, &
        &               dmin,ilim,wlim,halo)

    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell.
    ! DELX  grid-cell spacing array. LENGTH(DELX) == +1 if 
    !       spacing is uniform .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! FHAT  grid-cell re-con. array. FHAT is an array with
    !       SIZE = MDOF-by-NVAR-by-NPOS-1 .
    ! EDGE  edge-centred interp. for function-value. EDGE
    !       is an array with SIZE = NVAR-by-NPOS .
    ! OSCL  grid-cell oscil. dof.'s. OSCL is an array with
    !       SIZE = +2  -by-NVAR-by-NPOS-1 .
    ! DMIN  min. grid-cell spacing thresh . 
    ! ILIM  cell slope-limiting selection .
    ! WLIM  wall slope-limiting selection .
    ! HALO  width of re-con. stencil, symmetric about mid. .
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)  :: npos,nvar,ndof
        real*8 , intent(in)  :: dmin
        real*8 , intent(out) :: fhat(:,:,:)
        real*8 , intent(in)  :: oscl(:,:,:)
        real*8 , intent(in)  :: delx(:)
        real*8 , intent(in)  :: fdat(:,:,:)
        real*8 , intent(in)  :: edge(:,:)
        integer, intent(in)  :: ilim,wlim,halo

    !------------------------------------------- variables !
        integer :: ipos,ivar,iill,iirr,head,tail
        real*8  :: ff00,ffll,ffrr,hh00,hhll,hhrr
        integer :: mono
        real*8  :: fell,ferr
        real*8  :: dfds(-1:+1)
        real*8  :: wval(+1:+2)
        real*8  :: uhat(+1:+3)
        real*8  :: lhat(+1:+3)
        
        head = +1; tail = npos - 1

        if (npos.eq.2) then
    !----- default to reduced order if insufficient points !
        do  ivar = +1, nvar
            fhat(1,ivar,+1) = &
        &   fdat(1,ivar,+1)
            fhat(2,ivar,+1) = 0.d0
            fhat(3,ivar,+1) = 0.d0
        end do
        end if

        if (npos.le.2) return

    !------------------- reconstruct function on each cell !

        uhat = +0.d+0
        lhat = +0.d+0

        do  ipos = +1 , npos-1

        iill = max(head,ipos-1)
        iirr = min(tail,ipos+1)

        do  ivar = +1 , nvar-0

    !----------------------------- cell mean + edge values !
    
            ff00 = fdat(1,ivar,ipos)
            ffll = fdat(1,ivar,iill)
            ffrr = fdat(1,ivar,iirr)

            fell = edge(ivar,ipos+0)
            ferr = edge(ivar,ipos+1)

    !----------------------------- calc. LL/00/RR gradient !
    
            if (size(delx).gt.+1) then

            hh00 = delx(ipos)
            hhll = delx(iill)
            hhrr = delx(iirr)

            call plsv (ffll,hhll,ff00, &
    &                  hh00,ffrr,hhrr, &
    &                  dfds)
            else
            
            call plsc (ffll,ff00,ffrr, &
    &                  dfds)
    
            end if

    !----------------------------- calc. cell-wise profile !
    
            select case(ilim)
            case (null_limit)

    !----------------------------- calc. unlimited profile !   
               
            call ppmfn(ff00,ffll,ffrr, &
    &                  fell,ferr,dfds, &
    &                  uhat,lhat,mono)
            
    !----------------------------- pref. unlimited profile !

            wval(1) = +1.d+0
            wval(2) = +0.d+0
            
            case (mono_limit)
            
    !----------------------------- calc. monotonic profile ! 
             
            call ppmfn(ff00,ffll,ffrr, &
    &                  fell,ferr,dfds, &
    &                  uhat,lhat,mono)
            
    !----------------------------- pref. monotonic profile !

            wval(1) = +0.d+0
            wval(2) = +1.d+0
            
            case (weno_limit)

    !----------------------------- calc. unlimited profile !  
            
            call ppmfn(ff00,ffll,ffrr, &
    &                  fell,ferr,dfds, &
    &                  uhat,lhat,mono)

            if (mono.gt.+0) then
  
    !----------------------------- calc. WENO-type weights ! 
     
            call wenoi(npos,delx,oscl, &
    &                  ipos,ivar,halo, &
    &                  wlim,wval)
            
            else
                
    !----------------------------- pref. unlimited profile !

            wval(1) = +1.d+0
            wval(2) = +0.d+0
            
            end if
            
            end select
 
    !----------------------------- blend "null" and "mono" !
               
            fhat(1,ivar,ipos) = &
    &       wval(1) * uhat(1) + &
    &       wval(2) * lhat(1)
            fhat(2,ivar,ipos) = &
    &       wval(1) * uhat(2) + &
    &       wval(2) * lhat(2)
            fhat(3,ivar,ipos) = &
    &       wval(1) * uhat(3) + &
    &       wval(2) * lhat(3)

        end do

        end do
        
        return

    end  subroutine

    !--------- assemble piecewise parabolic reconstruction !
    
    pure subroutine ppmfn(ff00,ffll,ffrr,fell,&
        &                 ferr,dfds,uhat,lhat,&
        &                 mono)

    !
    ! FF00  centred grid-cell mean.
    ! FFLL  left -biased grid-cell mean.
    ! FFRR  right-biased grid-cell mean.
    ! FELL  left -biased edge interp.
    ! FERR  right-biased edge interp.
    ! DFDS  piecewise linear gradients in local co-ord.'s.
    !       DFDS(+0) is a centred, slope-limited estimate,
    !       DFDS(-1), DFDS(+1) are left- and right-biased
    !       estimates (unlimited). 
    ! UHAT  unlimited PPM reconstruction coefficients .
    ! LHAT  monotonic PPM reconstruction coefficients .
    ! MONO  slope-limiting indicator, MONO > +0 if some 
    !       limiting has occured .
    !

        implicit none

    !------------------------------------------- arguments !
        real*8 , intent(in)    :: ff00
        real*8 , intent(in)    :: ffll,ffrr
        real*8 , intent(inout) :: fell,ferr
        real*8 , intent(in)    :: dfds(-1:+1)
        real*8 , intent(out)   :: uhat(+1:+3)
        real*8 , intent(out)   :: lhat(+1:+3)
        integer, intent(out)   :: mono
          
    !------------------------------------------- variables !
        real*8  :: turn
        
        mono  = 0
        
    !-------------------------------- "null" slope-limiter !

        uhat( 1 ) = &
    & + (3.0d+0 / 2.0d+0) * ff00 &
    & - (1.0d+0 / 4.0d+0) *(ferr+fell)
        uhat( 2 ) = &
    & + (1.0d+0 / 2.0d+0) *(ferr-fell)
        uhat( 3 ) = &
    & - (3.0d+0 / 2.0d+0) * ff00 &
    & + (3.0d+0 / 4.0d+0) *(ferr+fell)
    
    !-------------------------------- "mono" slope-limiter !
        
        if((ffrr - ff00) * & 
    &      (ff00 - ffll) .lt. 0.d+0) then

    !----------------------------------- "flatten" extrema !
    
            mono = +1

            lhat(1) = ff00
            lhat(2) = 0.d0
            lhat(3) = 0.d0
             
            return
              
        end if

    !----------------------------------- limit edge values !
    
        if((ffll - fell) * &
    &      (fell - ff00) .le. 0.d+0) then

            mono = +1
       
            fell = ff00 - dfds(0)

        end if

        if((ffrr - ferr) * &
    &      (ferr - ff00) .le. 0.d+0) then

            mono = +1

            ferr = ff00 + dfds(0)
         
        end if
   
    !----------------------------------- update ppm coeff. !

        lhat( 1 ) = &
    & + (3.0d+0 / 2.0d+0) * ff00 &
    & - (1.0d+0 / 4.0d+0) *(ferr+fell)
        lhat( 2 ) = &
    & + (1.0d+0 / 2.0d+0) *(ferr-fell)
        lhat( 3 ) = &
    & - (3.0d+0 / 2.0d+0) * ff00 &
    & + (3.0d+0 / 4.0d+0) *(ferr+fell)

    !----------------------------------- limit cell values !
    
        if (abs(lhat(3)) .gt. & 
    &       abs(lhat(2))*.5d+0) then

        turn = -0.5d+0 * lhat(2) &
    &                  / lhat(3)

        if ((turn .ge. -1.d+0)&
    &  .and.(turn .le. +0.d+0)) then

        mono =   +2

    !--------------------------- push TURN onto lower edge !

        ferr =   +3.0d+0  * ff00 &
    &            -2.0d+0  * fell

        lhat( 1 ) = &
    & + (3.0d+0 / 2.0d+0) * ff00 &
    & - (1.0d+0 / 4.0d+0) *(ferr+fell)
        lhat( 2 ) = &
    & + (1.0d+0 / 2.0d+0) *(ferr-fell)
        lhat( 3 ) = &
    & - (3.0d+0 / 2.0d+0) * ff00 &
    & + (3.0d+0 / 4.0d+0) *(ferr+fell)

        else &
    &   if ((turn .gt. +0.d+0)&
    &  .and.(turn .le. +1.d+0)) then

        mono =   +2

    !--------------------------- push TURN onto upper edge !
    
        fell =   +3.0d+0  * ff00 &
    &            -2.0d+0  * ferr

        lhat( 1 ) = &
    & + (3.0d+0 / 2.0d+0) * ff00 &
    & - (1.0d+0 / 4.0d+0) *(ferr+fell)
        lhat( 2 ) = &
    & + (1.0d+0 / 2.0d+0) *(ferr-fell)
        lhat( 3 ) = &
    & - (3.0d+0 / 2.0d+0) * ff00 &
    & + (3.0d+0 / 4.0d+0) *(ferr+fell)
   
        end if
      
        end if
   
        return
    
    end subroutine



