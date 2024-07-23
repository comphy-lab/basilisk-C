
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
    ! PBC.f90: setup polynomial B.C.'s at domain endpoints.
    !
    ! Darren Engwirda 
    ! 09-Sep-2016
    ! de2363 [at] columbia [dot] edu
    !
    !

    subroutine pbc(npos,nvar,ndof,delx, &
        &          fdat,bcon,edge,dfdx, &
        &          iend,dmin)

    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell.
    ! DELX  grid-cell spacing array. LENGTH(DELX) == +1 if 
    !       spacing is uniform .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! BCON  boundary condition data for endpoint .
    ! EDGE  edge-centred interp. for function-value. EDGE
    !       is an array with SIZE = NVAR-by-NPOS .
    ! DFDX  edge-centred interp. for 1st-derivative. DFDX
    !       is an array with SIZE = NVAR-by-NPOS .
    ! IEND  domain endpoint, IEND < +0 for lower end-point 
    !       and IEND > +0 for upper endpoint .
    ! DMIN  min. grid-cell spacing thresh . 
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        real*8 , intent( in) :: delx(:)
        real*8 , intent( in) :: fdat(:,:,:)
        real*8 , intent(out) :: edge(:,:)
        real*8 , intent(out) :: dfdx(:,:)
        integer, intent( in) :: iend
        real*8 , intent( in) :: dmin
        type(rcon_ends), intent(in) :: bcon(:)

    !------------------------------------------- variables !
        integer :: ivar,nlse,nval,nslp

        nlse = 0 ; nval = 0 ; nslp = 0

        do  ivar = +1, nvar

            select case (bcon(ivar)%bcopt)
    !------------------------------------------- find BC's !
            case(bcon_loose)
                nlse = nlse + 1
            
            case(bcon_value)
                nval = nval + 1
            
            case(bcon_slope)
                nslp = nslp + 1

            end select

        end do

    !---------------------------- setup "lower" conditions !

        if (iend.lt.+0) then

        if (nlse.gt.+0) then
    !---------------------------- setup "unset" conditions !
        call lbc(npos,nvar,ndof, &
    &            delx,fdat,bcon, &
    &            bcon_loose    , &
    &            edge,dfdx,dmin)

        end if

        if (nval.gt.+0) then
    !---------------------------- setup "value" conditions !
        call lbc(npos,nvar,ndof, &
    &            delx,fdat,bcon, &
    &            bcon_value    , &
    &            edge,dfdx,dmin)

        end if

        if (nslp.gt.+0) then
    !---------------------------- setup "slope" conditions !
        call lbc(npos,nvar,ndof, &
    &            delx,fdat,bcon, &
    &            bcon_slope    , &
    &            edge,dfdx,dmin)

        end if
    
        end if
        
    !---------------------------- setup "upper" conditions !

        if (iend.gt.+0) then

        if (nlse.gt.+0) then
    !---------------------------- setup "unset" conditions !
        call ubc(npos,nvar,ndof, &
    &            delx,fdat,bcon, &
    &            bcon_loose    , &
    &            edge,dfdx,dmin)

        end if

        if (nval.gt.+0) then
    !---------------------------- setup "value" conditions !
        call ubc(npos,nvar,ndof, &
    &            delx,fdat,bcon, &
    &            bcon_value    , &
    &            edge,dfdx,dmin)

        end if

        if (nslp.gt.+0) then
    !---------------------------- setup "slope" conditions !
        call ubc(npos,nvar,ndof, &
    &            delx,fdat,bcon, &
    &            bcon_slope    , &
    &            edge,dfdx,dmin)

        end if
    
        end if
        
        return

    end  subroutine
       
    ! LBC: impose a single B.C.-type at the lower endpoint !
    
    subroutine lbc(npos,nvar,ndof,delx, &
        &          fdat,bcon,bopt,edge, &
        &          dfdx,dmin)

    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell.
    ! DELX  grid-cell spacing array. LENGTH(DELX) == +1 if 
    !       spacing is uniform .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! BCON  boundary condition data for endpoint .
    ! EDGE  edge-centred interp. for function-value. EDGE
    !       is an array with SIZE = NVAR-by-NPOS .
    ! DFDX  edge-centred interp. for 1st-derivative. DFDX
    !       is an array with SIZE = NVAR-by-NPOS .
    ! DMIN  min. grid-cell spacing thresh . 
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        integer, intent( in) :: bopt
        real*8 , intent( in) :: delx(:)
        real*8 , intent( in) :: fdat(:,:,:)
        real*8 , intent(out) :: edge(:,:)
        real*8 , intent(out) :: dfdx(:,:)
        real*8 , intent( in) :: dmin
        type(rcon_ends), intent(in) :: bcon(:)

    !------------------------------------------- variables !
        integer :: ivar,idof,isel, &
        &          head,tail,nsel
        logical :: okay        
        real*8  :: xhat
        real*8  :: delh(-1:+1)
        real*8  :: xmap(-1:+2)
        real*8  :: bvec(+3,-1:+2)
        real*8  :: gvec(+3,-1:+2)
        real*8  :: cmat(+3,+3)
        real*8  :: fhat(+3, nvar)
        real*8  :: eval(-1:+2)
        real*8  :: gval(-1:+2)
        
        integer, parameter :: NSIZ = +3
        real*8 , parameter :: ZERO = +1.d-14  

        head = +2; tail = npos - 2

        if (size(delx).gt.+1) then

    !------------------ mean grid spacing about ii-th cell !

        xhat = max(delx(head),dmin) * 0.5d+0

    !------------------ grid spacing for all stencil cells !

        delh(-1) = delx(head-1)
        delh(+0) = delx(head+0)
        delh(+1) = delx(head+1)

        else
        
    !------------------ mean grid spacing about ii-th cell !

        xhat = max(delx(  +1),dmin) * 0.5d+0

    !------------------ grid spacing for all stencil cells !

        delh(-1) = delx(    +1)
        delh(+0) = delx(    +1)
        delh(+1) = delx(    +1)
        
        end if

    !---------- local coordinate mapping for stencil edges !

        xmap(-1) =-(delh(-1) + &
        &           delh(+0)*0.5d0)/xhat
        xmap(+0) = -1.d0
        xmap(+1) = +1.d0
        xmap(+2) = (delh(+1) + &
        &           delh(+0)*0.5d0)/xhat

    !------------ linear system: lhs reconstruction matrix !

        select case(bopt )
        case( bcon_loose )

        call bfun1d(-1,+3,xmap(-1),bvec(:,-1))
        call bfun1d(-1,+3,xmap(+0),bvec(:,+0))
        call bfun1d(-1,+3,xmap(+1),bvec(:,+1))
        call bfun1d(-1,+3,xmap(+2),bvec(:,+2))

        do  idof = +1 , +3

            cmat(1,idof) = bvec(idof,+0) &
        &                - bvec(idof,-1)
            cmat(2,idof) = bvec(idof,+1) &
        &                - bvec(idof,+0)
            cmat(3,idof) = bvec(idof,+2) &
        &                - bvec(idof,+1)

        end do

        case( bcon_value )

        call bfun1d(+0,+3,xmap(-1),gvec(:,-1))
        
        call bfun1d(-1,+3,xmap(-1),bvec(:,-1))
        call bfun1d(-1,+3,xmap(+0),bvec(:,+0))
        call bfun1d(-1,+3,xmap(+1),bvec(:,+1))

        do  idof = +1 , +3

            cmat(1,idof) = bvec(idof,+0) &
        &                - bvec(idof,-1)
            cmat(2,idof) = bvec(idof,+1) &
        &                - bvec(idof,+0)

            cmat(3,idof) = gvec(idof,-1)

        end do

        case( bcon_slope )

        call bfun1d(+1,+3,xmap(-1),gvec(:,-1))
        
        call bfun1d(-1,+3,xmap(-1),bvec(:,-1))
        call bfun1d(-1,+3,xmap(+0),bvec(:,+0))
        call bfun1d(-1,+3,xmap(+1),bvec(:,+1))

        do  idof = +1 , +3

            cmat(1,idof) = bvec(idof,+0) &
        &                - bvec(idof,-1)
            cmat(2,idof) = bvec(idof,+1) &
        &                - bvec(idof,+0)

            cmat(3,idof) = gvec(idof,-1)

        end do

        end select

    !------------ linear system: rhs reconstruction vector !

        isel = 0 ; nsel = 0

        select case( bopt )
        case ( bcon_loose )

        do  ivar = +1, nvar

        if (bcon(ivar)%bcopt.eq.bcon_loose)  then

            isel = isel + 1
            nsel = nsel + 1

            fhat(1,isel) = delh(-1) * &
        &       fdat(1,ivar,head-1) / xhat
            fhat(2,isel) = delh(+0) * &
        &       fdat(1,ivar,head+0) / xhat
            fhat(3,isel) = delh(+1) * &
        &       fdat(1,ivar,head+1) / xhat
 
        end if
        
        end do       
 
        case ( bcon_value )
        
        do  ivar = +1, nvar

        if (bcon(ivar)%bcopt.eq.bcon_value)  then

            isel = isel + 1
            nsel = nsel + 1

            fhat(1,isel) = delh(-1) * &
        &       fdat(1,ivar,head-1) / xhat
            fhat(2,isel) = delh(+0) * &
        &       fdat(1,ivar,head+0) / xhat

            fhat(3,isel) = bcon(ivar)%value
 
        end if
        
        end do
        
        case ( bcon_slope )

        do  ivar = +1, nvar

        if (bcon(ivar)%bcopt.eq.bcon_slope)  then

            isel = isel + 1
            nsel = nsel + 1

            fhat(1,isel) = delh(-1) * &
        &       fdat(1,ivar,head-1) / xhat
            fhat(2,isel) = delh(+0) * &
        &       fdat(1,ivar,head+0) / xhat

            fhat(3,isel) = &
        &       bcon(ivar)%slope * xhat
 
        end if
        
        end do
        
        end select

    !------------------------- factor/solve linear systems !

        call slv_3x3(cmat,NSIZ,fhat , &
        &            NSIZ,nvar,       &
        &            ZERO*dmin,okay)

        if (okay .eqv..false.) then

#       ifdef __PPR_WARNMAT__
        
        write(*,*) &
        &   "WARNING::LBC-matrix-is-singular!"
        
#       endif

        end if

        if (okay .eqv. .true.) then

    !------------- extrapolate values/slopes at lower edge !

        isel  = +0

        call bfun1d(+0,+3,xmap(-1),bvec(:,-1))
        call bfun1d(+0,+3,xmap(+0),bvec(:,+0))
        call bfun1d(+0,+3,xmap(+1),bvec(:,+1))
        
        call bfun1d(+1,+3,xmap(-1),gvec(:,-1))
        call bfun1d(+1,+3,xmap(+0),gvec(:,+0))
        call bfun1d(+1,+3,xmap(+1),gvec(:,+1))

        do  ivar = +1, nvar

        if (bcon(ivar)%bcopt.eq.bopt)  then

            isel = isel  + 1

            eval(-1) = dot_product( &
        &       bvec(:,-1),fhat(:,isel))
            eval(+0) = dot_product( &
        &       bvec(:,+0),fhat(:,isel))        
            eval(+1) = dot_product( &
        &       bvec(:,+1),fhat(:,isel))
        
            gval(-1) = dot_product( &
        &       gvec(:,-1),fhat(:,isel))
            gval(+0) = dot_product( &
        &       gvec(:,+0),fhat(:,isel))        
            gval(+1) = dot_product( &
        &       gvec(:,+1),fhat(:,isel))

            edge(ivar,head-1) = eval(-1)
            edge(ivar,head+0) = eval(+0)
            edge(ivar,head+1) = eval(+1)

            dfdx(ivar,head-1) = gval(-1) &
        &                     / xhat
            dfdx(ivar,head+0) = gval(+0) &
        &                     / xhat
            dfdx(ivar,head+1) = gval(+1) &
        &                     / xhat

        end if

        end do
   
        else

    !------------- low-order if re-con. matrix is singular !

        do  ivar = +1, nvar

        if (bcon(ivar)%bcopt.eq.bopt)  then

            eval(-1) = &
        &   fdat(1,ivar,head-1) * 1.d0
            eval(+0) = &
        &   fdat(1,ivar,head-1) * .5d0 + &
        &   fdat(1,ivar,head+0) * .5d0
            eval(+1) = &
        &   fdat(1,ivar,head+0) * .5d0 + &
        &   fdat(1,ivar,head+1) * .5d0
   
            gval(-1) = &
        &   fdat(1,ivar,head+0) * .5d0 - &
        &   fdat(1,ivar,head-1) * .5d0
            gval(+0) = &
        &   fdat(1,ivar,head+0) * .5d0 - &
        &   fdat(1,ivar,head-1) * .5d0
            gval(+1) = &
        &   fdat(1,ivar,head+1) * .5d0 - &
        &   fdat(1,ivar,head+0) * .5d0
   
            edge(ivar,head-1) = eval(-1)
            edge(ivar,head+0) = eval(+0)
            edge(ivar,head+1) = eval(+1)

            dfdx(ivar,head-1) = gval(-1) &
        &                     / xhat
            dfdx(ivar,head+0) = gval(+0) &
        &                     / xhat
            dfdx(ivar,head+1) = gval(+1) &
        &                     / xhat

        end if
        
        end do
    
        end if
        
        return

    end  subroutine
    
    ! UBC: impose a single B.C.-type at the upper endpoint !
    
    subroutine ubc(npos,nvar,ndof,delx, &
        &          fdat,bcon,bopt,edge, &
        &          dfdx,dmin)

    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell.
    ! DELX  grid-cell spacing array. LENGTH(DELX) == +1 if 
    !       spacing is uniform .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! BCON  boundary condition data for endpoint .
    ! EDGE  edge-centred interp. for function-value. EDGE
    !       is an array with SIZE = NVAR-by-NPOS .
    ! DFDX  edge-centred interp. for 1st-derivative. DFDX
    !       is an array with SIZE = NVAR-by-NPOS .
    ! DMIN  min. grid-cell spacing thresh . 
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        integer, intent( in) :: bopt
        real*8 , intent( in) :: delx(:)
        real*8 , intent( in) :: fdat(:,:,:)
        real*8 , intent(out) :: edge(:,:)
        real*8 , intent(out) :: dfdx(:,:)
        real*8 , intent( in) :: dmin
        type(rcon_ends), intent(in) :: bcon(:)

    !------------------------------------------- variables !
        integer :: ivar,idof,isel, &
        &          head,tail,nsel
        logical :: okay
        real*8  :: xhat
        real*8  :: delh(-1:+1)
        real*8  :: xmap(-1:+2)
        real*8  :: bvec(+3,-1:+2)
        real*8  :: gvec(+3,-1:+2)
        real*8  :: cmat(+3,+3)
        real*8  :: fhat(+3, nvar)
        real*8  :: eval(-1:+2)
        real*8  :: gval(-1:+2)

        integer, parameter :: NSIZ = +3
        real*8 , parameter :: ZERO = +1.d-14

        head = +2; tail = npos - 2

        if (size(delx).gt.+1) then

    !------------------ mean grid spacing about ii-th cell !

        xhat = max(delx(tail),dmin) * 0.5d+0

    !------------------ grid spacing for all stencil cells !

        delh(-1) = delx(tail-1)
        delh(+0) = delx(tail+0)
        delh(+1) = delx(tail+1)

        else
        
    !------------------ mean grid spacing about ii-th cell !

        xhat = max(delx(  +1),dmin) * 0.5d+0

    !------------------ grid spacing for all stencil cells !

        delh(-1) = delx(    +1)
        delh(+0) = delx(    +1)
        delh(+1) = delx(    +1)
        
        end if

    !---------- local coordinate mapping for stencil edges !

        xmap(-1) =-(delh(-1) + &
        &           delh(+0)*0.5d0)/xhat
        xmap(+0) = -1.d0
        xmap(+1) = +1.d0
        xmap(+2) = (delh(+1) + &
        &           delh(+0)*0.5d0)/xhat

    !------------ linear system: lhs reconstruction matrix !

        select case(bopt )
        case( bcon_loose )

        call bfun1d(-1,+3,xmap(-1),bvec(:,-1))
        call bfun1d(-1,+3,xmap(+0),bvec(:,+0))
        call bfun1d(-1,+3,xmap(+1),bvec(:,+1))
        call bfun1d(-1,+3,xmap(+2),bvec(:,+2))

        do  idof = +1 , +3

            cmat(1,idof) = bvec(idof,+0) &
        &                - bvec(idof,-1)
            cmat(2,idof) = bvec(idof,+1) &
        &                - bvec(idof,+0)
            cmat(3,idof) = bvec(idof,+2) &
        &                - bvec(idof,+1)

        end do

        case( bcon_value )

        call bfun1d(+0,+3,xmap(+2),gvec(:,+2))
        
        call bfun1d(-1,+3,xmap(+0),bvec(:,+0))
        call bfun1d(-1,+3,xmap(+1),bvec(:,+1))
        call bfun1d(-1,+3,xmap(+2),bvec(:,+2))

        do  idof = +1 , +3

            cmat(1,idof) = bvec(idof,+1) &
        &                - bvec(idof,+0)
            cmat(2,idof) = bvec(idof,+2) &
        &                - bvec(idof,+1)

            cmat(3,idof) = gvec(idof,+2)

        end do

        case( bcon_slope )

        call bfun1d(+1,+3,xmap(+2),gvec(:,+2))
        
        call bfun1d(-1,+3,xmap(+0),bvec(:,+0))
        call bfun1d(-1,+3,xmap(+1),bvec(:,+1))
        call bfun1d(-1,+3,xmap(+2),bvec(:,+2))

        do  idof = +1 , +3

            cmat(1,idof) = bvec(idof,+1) &
        &                - bvec(idof,+0)
            cmat(2,idof) = bvec(idof,+2) &
        &                - bvec(idof,+1)

            cmat(3,idof) = gvec(idof,+2)

        end do

        end select

    !------------ linear system: rhs reconstruction vector !

        isel = 0 ; nsel = 0

        select case( bopt )
        case ( bcon_loose )

        do  ivar = +1, nvar

        if (bcon(ivar)%bcopt.eq.bcon_loose)  then

            isel = isel + 1
            nsel = nsel + 1

            fhat(1,isel) = delh(-1) * &
        &       fdat(1,ivar,tail-1) / xhat
            fhat(2,isel) = delh(+0) * &
        &       fdat(1,ivar,tail+0) / xhat
            fhat(3,isel) = delh(+1) * &
        &       fdat(1,ivar,tail+1) / xhat
 
        end if
        
        end do       
 
        case ( bcon_value )
        
        do  ivar = +1, nvar

        if (bcon(ivar)%bcopt.eq.bcon_value)  then

            isel = isel + 1
            nsel = nsel + 1

            fhat(1,isel) = delh(+0) * &
        &       fdat(1,ivar,tail+0) / xhat
            fhat(2,isel) = delh(+1) * &
        &       fdat(1,ivar,tail+1) / xhat

            fhat(3,isel) = bcon(ivar)%value
 
        end if
        
        end do
        
        case ( bcon_slope )

        do  ivar = +1, nvar

        if (bcon(ivar)%bcopt.eq.bcon_slope)  then

            isel = isel + 1
            nsel = nsel + 1

            fhat(1,isel) = delh(+0) * &
        &       fdat(1,ivar,tail+0) / xhat
            fhat(2,isel) = delh(+1) * &
        &       fdat(1,ivar,tail+1) / xhat

            fhat(3,isel) = &
        &       bcon(ivar)%slope * xhat
 
        end if
        
        end do
        
        end select

    !------------------------- factor/solve linear systems !

        call slv_3x3(cmat,NSIZ,fhat , &
        &            NSIZ,nvar,       &
        &            ZERO*dmin,okay)

        if (okay .eqv..false.) then

#       ifdef __PPR_WARNMAT__
        
        write(*,*) &
        &   "WARNING::UBC-matrix-is-singular!"
        
#       endif

        end if

        if (okay .eqv. .true.) then

    !------------- extrapolate values/slopes at lower edge !

        isel  = +0

        call bfun1d(+0,+3,xmap(+0),bvec(:,+0))
        call bfun1d(+0,+3,xmap(+1),bvec(:,+1))
        call bfun1d(+0,+3,xmap(+2),bvec(:,+2))
        
        call bfun1d(+1,+3,xmap(+0),gvec(:,+0))
        call bfun1d(+1,+3,xmap(+1),gvec(:,+1))
        call bfun1d(+1,+3,xmap(+2),gvec(:,+2))

        do  ivar = +1, nvar

        if (bcon(ivar)%bcopt.eq.bopt)  then

            isel = isel  + 1

            eval(+0) = dot_product( &
        &       bvec(:,+0),fhat(:,isel))
            eval(+1) = dot_product( &
        &       bvec(:,+1),fhat(:,isel))        
            eval(+2) = dot_product( &
        &       bvec(:,+2),fhat(:,isel))
        
            gval(+0) = dot_product( &
        &       gvec(:,+0),fhat(:,isel))
            gval(+1) = dot_product( &
        &       gvec(:,+1),fhat(:,isel))        
            gval(+2) = dot_product( &
        &       gvec(:,+2),fhat(:,isel))

            edge(ivar,tail+0) = eval(+0)
            edge(ivar,tail+1) = eval(+1)
            edge(ivar,tail+2) = eval(+2)

            dfdx(ivar,tail+0) = gval(+0) &
        &                     / xhat
            dfdx(ivar,tail+1) = gval(+1) &
        &                     / xhat
            dfdx(ivar,tail+2) = gval(+2) &
        &                     / xhat

        end if

        end do
   
        else

    !------------- low-order if re-con. matrix is singular !

        do  ivar = +1, nvar

        if (bcon(ivar)%bcopt.eq.bopt)  then

            eval(+0) = &
        &   fdat(1,ivar,tail-1) * .5d0 + &
        &   fdat(1,ivar,tail+0) * .5d0
            eval(+1) = &
        &   fdat(1,ivar,tail+0) * .5d0 + &
        &   fdat(1,ivar,tail+1) * .5d0
            eval(+2) = &
        &   fdat(1,ivar,tail+1) * 1.d0
   
            gval(+0) = &
        &   fdat(1,ivar,tail+0) * .5d0 - &
        &   fdat(1,ivar,tail-1) * .5d0
            gval(+1) = &
        &   fdat(1,ivar,tail+1) * .5d0 - &
        &   fdat(1,ivar,tail+0) * .5d0
            gval(+2) = &
        &   fdat(1,ivar,tail+1) * .5d0 - &
        &   fdat(1,ivar,tail+0) * .5d0
   
            edge(ivar,tail+0) = eval(+0)
            edge(ivar,tail+1) = eval(+1)
            edge(ivar,tail+2) = eval(+2)

            dfdx(ivar,tail+0) = gval(+0) &
        &                     / xhat
            dfdx(ivar,tail+1) = gval(+1) &
        &                     / xhat
            dfdx(ivar,tail+2) = gval(+2) &
        &                     / xhat

        end if
        
        end do
    
        end if
        
        return

    end  subroutine
    
    
    
