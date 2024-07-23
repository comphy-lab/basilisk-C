    
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
    ! OSCL1D.f90: "oscillation-indicators" for WENO interp.
    !
    ! Darren Engwirda 
    ! 08-Sep-2016
    ! de2363 [at] columbia [dot] edu
    !
    !
    
    pure subroutine oscli (npos,nvar,ndof,delx,&
        &                  fdat,oscl,dmin)

    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell .
    ! DELX  (constant) grid-cell spacing. LENGTH(DELX)==+1 .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! OSCL  grid-cell oscil. dof.'s. OSCL is an array with
    !       SIZE = +2  -by-NVAR-by-NPOS-1 .
    ! DMIN  min. grid-cell spacing thresh .
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        real*8 , intent( in) :: dmin
        real*8 , intent( in) :: delx(:)
        real*8 , intent( in) :: fdat(:,:,:)
        real*8 , intent(out) :: oscl(:,:,:)

    !------------------------------------------- variables !
        integer :: ivar,ipos
        
        if (npos.lt.3) then
    !------------------------------- at least 3 grid-cells !
        do  ipos = +1, npos-1
        do  ivar = +1, nvar-0
            oscl(1,ivar,ipos) = +0.d0
            oscl(2,ivar,ipos) = +0.d0
        end do
        end do
        end if
        
        if (npos.lt.3) return
        if (nvar.lt.1) return
        if (ndof.lt.1) return

        if (size(delx).gt.+1) then

    !------------------------------- variable grid-spacing !

            call osclv(npos,nvar,ndof,delx, &
        &              fdat,oscl,dmin)
        
        else

    !------------------------------- constant grid-spacing !
        
            call osclc(npos,nvar,ndof,delx, &
        &              fdat,oscl,dmin)
                
        end if

        return
        
    end  subroutine
    
    pure subroutine osclv (npos,nvar,ndof,delx,&
        &                  fdat,oscl,dmin)

    !
    ! *this is the variable grid-spacing variant .
    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell .
    ! DELX  (variable) grid-cell spacing. LENGTH(DELX)!=+1 .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! OSCL  grid-cell oscil. dof.'s. OSCL is an array with
    !       SIZE = +2  -by-NVAR-by-NPOS-1 .
    ! DMIN  min. grid-cell spacing thresh .
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        real*8 , intent( in) :: dmin
        real*8 , intent( in) :: delx(:)
        real*8 , intent( in) :: fdat(:,:,:)
        real*8 , intent(out) :: oscl(:,:,:)

    !------------------------------------------- variables !
        integer :: head,tail
        integer :: ipos,ivar
        real*8  :: hhll,hhcc,hhrr
        real*8  :: hhmm,hhrc,hhlc
        real*8  :: cmat(2,3)

        head = +1 ; tail = npos-1

    !--------------------------------------- centred point !

        do  ipos = head+1, tail-1

        hhll = max(delx(ipos-1),dmin)
        hhcc = max(delx(ipos+0),dmin)
        hhrr = max(delx(ipos+1),dmin)

        hhrc = hhrr + hhcc
        hhlc = hhll + hhcc
        hhmm = hhll + hhcc + hhrr

        cmat(1,1) = -(hhcc+2.d0*hhrr)/(hhlc*hhmm)
        cmat(1,2) = -(hhll-hhrr)* &
        &          (3.d0*hhcc+2.d0*(hhll+hhrr))/&
        &            (hhlc*hhrc*hhmm)
        cmat(1,3) = +(hhcc+2.d0*hhll)/(hhrc*hhmm)

        cmat(2,1) = +3.d0/(hhlc*hhmm)
        cmat(2,2) = -3.d0*(2.d0*hhcc+hhll+hhrr)/&
        &            (hhlc*hhrc*hhmm)
        cmat(2,3) = +3.d0/(hhrc*hhmm)

        do  ivar = 1, nvar

            oscl(1,ivar,ipos) = +1.d0 * ( &
        & + cmat(1,1)*fdat(1,ivar,ipos-1) &
        & + cmat(1,2)*fdat(1,ivar,ipos+0) &
        & + cmat(1,3)*fdat(1,ivar,ipos+1) )

            oscl(2,ivar,ipos) = +2.d0 * ( &
        & + cmat(2,1)*fdat(1,ivar,ipos-1) &
        & + cmat(2,2)*fdat(1,ivar,ipos+0) &
        & + cmat(2,3)*fdat(1,ivar,ipos+1) )

        end do
            
        end do

    !-------------------------------------- lower endpoint !

        hhll = max(delx(head+0),dmin)
        hhcc = max(delx(head+1),dmin)
        hhrr = max(delx(head+2),dmin)

        cmat(1,1) = -2.d0 / (hhll+hhcc)
        cmat(1,2) = +2.d0 / (hhll+hhcc)

        do  ivar = 1, nvar

            oscl(1,ivar,head) = &
        & + cmat(1,1)*fdat(1,ivar,head+0) &
        & + cmat(1,2)*fdat(1,ivar,head+1)
             
            oscl(2,ivar,head) = +0.d0

        end do

    !-------------------------------------- upper endpoint !

        hhll = max(delx(tail-2),dmin)
        hhcc = max(delx(tail-1),dmin)
        hhrr = max(delx(tail-0),dmin)

        cmat(1,2) = -2.d0 / (hhrr+hhcc)
        cmat(1,3) = +2.d0 / (hhrr+hhcc)

        do  ivar = 1, nvar

            oscl(1,ivar,tail) = &
        & + cmat(1,2)*fdat(1,ivar,tail-1) &
        & + cmat(1,3)*fdat(1,ivar,tail+0)

            oscl(2,ivar,tail) = +0.d0

        end do
        
        return

    end  subroutine
    
    pure subroutine osclc (npos,nvar,ndof,delx,&
        &                  fdat,oscl,dmin)

    !
    ! *this is the constant grid-spacing variant .
    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell .
    ! DELX  (constant) grid-cell spacing. LENGTH(DELX)==+1 .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! OSCL  grid-cell oscil. dof.'s. OSCL is an array with
    !       SIZE = +2  -by-NVAR-by-NPOS-1 .
    ! DMIN  min. grid-cell spacing thresh .
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        real*8 , intent( in) :: dmin
        real*8 , intent( in) :: delx(1)
        real*8 , intent( in) :: fdat(:,:,:)
        real*8 , intent(out) :: oscl(:,:,:)

    !------------------------------------------- variables !
        integer :: head,tail,ipos,ivar
        
        head = +1; tail = npos - 1

    !-------------------------------------- centred points !

        do  ipos = 2, npos-2
        do  ivar = 1, nvar-0

            oscl(1,ivar,ipos) = &
        & + .25d+0 * fdat(1,ivar,ipos+1) &
        & - .25d+0 * fdat(1,ivar,ipos-1)
        
            oscl(2,ivar,ipos) = &
        & + .25d+0 * fdat(1,ivar,ipos+1) &
        & - .50d+0 * fdat(1,ivar,ipos+0) &
        & + .25d+0 * fdat(1,ivar,ipos-1)

        end do
        end do

    !-------------------------------------- lower endpoint !

        do  ivar = 1, nvar

            oscl(1,ivar,head) = &
        & + .50d+0 * fdat(1,ivar,head+1) &
        & - .50d+0 * fdat(1,ivar,head+0)
        
            oscl(2,ivar,head) = +0.d0

        end do

    !-------------------------------------- upper endpoint !

        do  ivar = 1, nvar

            oscl(1,ivar,tail) = &
        & + .50d+0 * fdat(1,ivar,tail+0) &
        & - .50d+0 * fdat(1,ivar,tail-1)
            
            oscl(2,ivar,tail) = +0.d0

        end do
        
        return

    end  subroutine
    
    

