
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
    ! P3E.f90: set edge estimates via degree-3 polynomials.
    !
    ! Darren Engwirda 
    ! 09-Sep-2016
    ! de2363 [at] columbia [dot] edu
    !
    !
    
    subroutine p3e(npos,nvar,ndof,delx, &
        &          fdat,bclo,bchi,edge, &
        &          dfdx,opts,dmin)

    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell.
    ! DELX  grid-cell spacing array. LENGTH(DELX) == +1 if 
    !       spacing is uniform .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! BCLO  boundary condition at lower endpoint.
    ! BCHI  boundary condition at upper endpoint.
    ! EDGE  edge-centred interp. for function-value. EDGE
    !       is an array with SIZE = NVAR-by-NPOS .
    ! DFDX  edge-centred interp. for 1st-derivative. DFDX
    !       is an array with SIZE = NVAR-by-NPOS .
    ! OPTS  method parameters. See RCON-OPTS for details .
    ! DMIN  min. grid-cell spacing thresh . 
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        real*8 , intent( in) :: delx(:)
        real*8 , intent( in) :: fdat(:,:,:)
        type (rcon_ends), intent(in) :: bclo(:)
        type (rcon_ends), intent(in) :: bchi(:)
        real*8 , intent(out) :: edge(:,:)
        real*8 , intent(out) :: dfdx(:,:)
        real*8 , intent( in) :: dmin
        class(rcon_opts), intent(in) :: opts

    !------------------------------------------- variables !
        integer :: ipos,ivar,idof,head,tail
        logical :: okay
        real*8  :: xhat,fEPS
        real*8  :: delh(-2:+1)
        real*8  :: xmap(-2:+2)
        real*8  :: fhat(+4, nvar)
        real*8  :: ivec(+4,-2:+2)
        real*8  :: cmat(+4,+4)
        
        integer, parameter :: NSIZ = +4
        real*8 , parameter :: ZERO = 1.d-14  
   
        head = +3 ; tail = npos - 2

        if (npos.le.4) then
    !----- default to reduced order if insufficient points !
            call p1e (npos,nvar,ndof, &
        &             delx,fdat,bclo, &
        &             bchi,edge,dfdx, &
        &             opts,dmin)
        end if

        if (npos.le.4) return

    !------ impose value/slope B.C.'s about lower endpoint !

        call pbc(npos,nvar,ndof,delx, &
        &        fdat,bclo,edge,dfdx, &
        &        -1  ,dmin)

    !------ impose value/slope B.C.'s about upper endpoint !

        call pbc(npos,nvar,ndof,delx, &
        &        fdat,bchi,edge,dfdx, &
        &        +1  ,dmin)

    ! Reconstruct edge-centred 4th-order polynomials. Com- !
    ! pute values/slopes at edges directly. Mid.-order ex- !
    ! trapolation at endpoints.                            !
    
        if (size(delx).eq.+1) then
        
            do  ipos = head , tail
    
    !--------------- reconstruction: constant grid-spacing !
            
            do  ivar = 1, nvar

                edge(ivar,ipos) = ( &
        &     -      1.d0 * &
        &       fdat(1,ivar,ipos-2) &
        &     +      7.d0 * &
        &       fdat(1,ivar,ipos-1) &
        &     +      7.d0 * &
        &       fdat(1,ivar,ipos+0) &
        &     -      1.d0 * &
        &       fdat(1,ivar,ipos+1) ) / 12.d0
            
                dfdx(ivar,ipos) = ( &
        &     +      1.d0 * &
        &       fdat(1,ivar,ipos-2) &
        &     -     15.d0 * &
        &       fdat(1,ivar,ipos-1) &
        &     +     15.d0 * &
        &       fdat(1,ivar,ipos+0) &
        &     -      1.d0 * &
        &       fdat(1,ivar,ipos+1) ) / 12.d0

                dfdx(ivar,ipos) = &
        &       dfdx(ivar,ipos) / delx(+1)

            end do
                      
            end do
        
        else
        
            fEPS     = ZERO * dmin

            do  ipos = head , tail
    
    !--------------- reconstruction: variable grid-spacing !
            
            delh(-2) = delx(ipos-2)
            delh(-1) = delx(ipos-1)
            delh(+0) = delx(ipos+0)
            delh(+1) = delx(ipos+1)

            xhat = .5d0 * max(delh(-1),dmin) + &
        &          .5d0 * max(delh(+0),dmin)

            xmap(-2) = -( delh(-2) &
        &              +  delh(-1) ) / xhat
            xmap(-1) = -  delh(-1)   / xhat
            xmap(+0) = +  0.d0
            xmap(+1) = +  delh(+0)   / xhat
            xmap(+2) = +( delh(+0) &
        &              +  delh(+1) ) / xhat
            
    !--------------------------- calc. integral basis vec. !

            do  idof = -2, +2

            ivec(1,idof) = &
        &       xmap(idof) ** 1 / 1.0d+0
            ivec(2,idof) = &
        &       xmap(idof) ** 2 / 2.0d+0
            ivec(3,idof) = &
        &       xmap(idof) ** 3 / 3.0d+0
            ivec(4,idof) = &
        &       xmap(idof) ** 4 / 4.0d+0
        
            end do

    !--------------------------- linear system: lhs matrix !

            do  idof = +1, +4

            cmat(1,idof) = ivec(idof,-1) &
        &                - ivec(idof,-2)
            cmat(2,idof) = ivec(idof,+0) &
        &                - ivec(idof,-1)
            cmat(3,idof) = ivec(idof,+1) &
        &                - ivec(idof,+0)
            cmat(4,idof) = ivec(idof,+2) &
        &                - ivec(idof,+1)

            end do

    !--------------------------- linear system: rhs vector !

            do  ivar = +1, nvar

            fhat(+1,ivar) = &
        &       delx(ipos-2) * &
        &   fdat(+1,ivar,ipos-2) / xhat
            fhat(+2,ivar) = &
        &       delx(ipos-1) * &
        &   fdat(+1,ivar,ipos-1) / xhat
            fhat(+3,ivar) = &
        &       delx(ipos+0) * &
        &   fdat(+1,ivar,ipos+0) / xhat
            fhat(+4,ivar) = &
        &       delx(ipos+1) * &
        &   fdat(+1,ivar,ipos+1) / xhat
        
            end do
            
    !------------------------- factor/solve linear systems !

            call slv_4x4(cmat,NSIZ,fhat, &
       &                 NSIZ,nvar,fEPS, &
       &                 okay)

            if (okay .eqv. .true.) then

            do  ivar =  +1, nvar

            edge(ivar,ipos) = fhat(1,ivar) 
                  
            dfdx(ivar,ipos) = fhat(2,ivar) &
        &                   / xhat
        
            end do

            else

    !------------------------- fallback if system singular !

#           ifdef __PPR_WARNMAT__
            
            write(*,*) &
    &   "WARNING::P3E - matrix-is-singular!"

#           endif

            do  ivar =  +1, nvar

            edge(ivar,ipos) =   &
        &   fdat(1,ivar,ipos-1) * 0.5d+0 + &
        &   fdat(1,ivar,ipos-0) * 0.5d+0
        
            dfdx(ivar,ipos) =   &
        &   fdat(1,ivar,ipos-0) * 1.0d+0 - &
        &   fdat(1,ivar,ipos-1) * 1.0d+0
                  
            dfdx(ivar,ipos) =   &
        &       dfdx(ivar,ipos) / xhat
        
            end do

            end if
        
            end do
        
        end if

        return
        
    end subroutine



