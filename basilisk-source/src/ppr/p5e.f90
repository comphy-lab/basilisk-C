
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
    ! P5E.f90: set edge estimates via degree-5 polynomials.
    !
    ! Darren Engwirda 
    ! 25-Mar-2019
    ! de2363 [at] columbia [dot] edu
    !
    !

    subroutine p5e(npos,nvar,ndof,delx, &
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
        real*8  :: delh(-3:+2)
        real*8  :: xmap(-3:+3)
        real*8  :: fhat(+6, nvar)
        real*8  :: ivec(+6,-3:+3)
        real*8  :: cmat(+6,+6)
        
        integer, parameter :: NSIZ = +6
        real*8 , parameter :: ZERO = 1.d-14   

        head = +4 ; tail = npos - 3

        if (npos.le.6) then
    !----- default to reduced order if insufficient points !
            call p3e (npos,nvar,ndof, &
        &             delx,fdat,bclo, &
        &             bchi,edge,dfdx, &
        &             opts,dmin)
        end if

        if (npos.le.6) return

    !------ impose value/slope B.C.'s about lower endpoint !

        call pbc(npos,nvar,ndof,delx, &
        &        fdat,bclo,edge,dfdx, &
        &        -1  ,dmin)

    !------ impose value/slope B.C.'s about upper endpoint !

        call pbc(npos,nvar,ndof,delx, &
        &        fdat,bchi,edge,dfdx, &
        &        +1  ,dmin)

    ! Reconstruct edge-centred 6th-order polynomials. Com- !
    ! pute values/slopes at edges directly. Mid.-order ex- !
    ! trapolation at endpoints.                            !
    
        if (size(delx).eq.+1) then
        
            do  ipos = head , tail
    
    !--------------- reconstruction: constant grid-spacing !
            
            do  ivar = 1, nvar

                edge(ivar,ipos) = &
        &     + ( 1.d0 / 60.d0) * & 
        &       fdat(1,ivar,ipos-3) &
        &     - ( 8.d0 / 60.d0) * &
        &       fdat(1,ivar,ipos-2) &
        &     + (37.d0 / 60.d0) * &
        &       fdat(1,ivar,ipos-1) &
        &     + (37.d0 / 60.d0) * &
        &       fdat(1,ivar,ipos+0) &
        &     - ( 8.d0 / 60.d0) * &
        &       fdat(1,ivar,ipos+1) &
        &     + ( 1.d0 / 60.d0) * &
        &       fdat(1,ivar,ipos+2)
            
                dfdx(ivar,ipos) = &
        &     - ( 1.d0 / 90.d0) * & 
        &       fdat(1,ivar,ipos-3) &
        &     + ( 5.d0 / 36.d0) * &
        &       fdat(1,ivar,ipos-2) &
        &     - (49.d0 / 36.d0) * &
        &       fdat(1,ivar,ipos-1) &
        &     + (49.d0 / 36.d0) * &
        &       fdat(1,ivar,ipos+0) &
        &     - ( 5.d0 / 36.d0) * &
        &       fdat(1,ivar,ipos+1) &
        &     + ( 1.d0 / 90.d0) * &
        &       fdat(1,ivar,ipos+2)

                dfdx(ivar,ipos) = &
                dfdx(ivar,ipos) / delx(+1)

            end do
                      
            end do
        
        else
        
            fEPS     = ZERO * dmin

            do  ipos = head , tail
    
    !--------------- reconstruction: variable grid-spacing !
            
            delh(-3) = &
        &       max(delx(ipos-3),dmin)
            delh(-2) = &
        &       max(delx(ipos-2),dmin)
            delh(-1) = &
        &       max(delx(ipos-1),dmin)
            delh(+0) = &
        &       max(delx(ipos+0),dmin)
            delh(+1) = &
        &       max(delx(ipos+1),dmin)
            delh(+2) = &
        &       max(delx(ipos+2),dmin)

            xhat = .5d0 * delh(-1) + &
        &          .5d0 * delh(+0)

            xmap(-3) = -( delh(-3) &
        &              +  delh(-2) &
        &              +  delh(-1) ) / xhat
            xmap(-2) = -( delh(-2) &
        &              +  delh(-1) ) / xhat
            xmap(-1) = -  delh(-1)   / xhat
            xmap(+0) = +  0.d0
            xmap(+1) = +  delh(+0)   / xhat
            xmap(+2) = +( delh(+0) &
        &              +  delh(+1) ) / xhat
            xmap(+3) = +( delh(+0) &
        &              +  delh(+1) &
        &              +  delh(+2) ) / xhat
            
    !--------------------------- calc. integral basis vec. !

            do  idof = -3, +3

            ivec(1,idof) = &
        &       xmap(idof) ** 1 / 1.0d+0
            ivec(2,idof) = &
        &       xmap(idof) ** 2 / 2.0d+0
            ivec(3,idof) = &
        &       xmap(idof) ** 3 / 3.0d+0
            ivec(4,idof) = &
        &       xmap(idof) ** 4 / 4.0d+0
            ivec(5,idof) = &
        &       xmap(idof) ** 5 / 5.0d+0
            ivec(6,idof) = &
        &       xmap(idof) ** 6 / 6.0d+0
        
            end do

    !--------------------------- linear system: lhs matrix !

            do  idof = +1, +6

            cmat(1,idof) = ivec(idof,-2) &
        &                - ivec(idof,-3)
            cmat(2,idof) = ivec(idof,-1) &
        &                - ivec(idof,-2)
            cmat(3,idof) = ivec(idof,+0) &
        &                - ivec(idof,-1)
            cmat(4,idof) = ivec(idof,+1) &
        &                - ivec(idof,+0)
            cmat(5,idof) = ivec(idof,+2) &
        &                - ivec(idof,+1)
            cmat(6,idof) = ivec(idof,+3) &
        &                - ivec(idof,+2)

            end do

    !--------------------------- linear system: rhs vector !

            do  ivar = +1, nvar

            fhat(+1,ivar) = &
        &       delx(ipos-3) * &
        &   fdat(+1,ivar,ipos-3) / xhat
            fhat(+2,ivar) = &
        &       delx(ipos-2) * &
        &   fdat(+1,ivar,ipos-2) / xhat
            fhat(+3,ivar) = &
        &       delx(ipos-1) * &
        &   fdat(+1,ivar,ipos-1) / xhat
            fhat(+4,ivar) = &
        &       delx(ipos+0) * &
        &   fdat(+1,ivar,ipos+0) / xhat
            fhat(+5,ivar) = &
        &       delx(ipos+1) * &
        &   fdat(+1,ivar,ipos+1) / xhat
            fhat(+6,ivar) = &
        &       delx(ipos+2) * &
        &   fdat(+1,ivar,ipos+2) / xhat
        
            end do
            
    !------------------------- factor/solve linear systems !

            call slv_6x6(cmat,NSIZ,fhat, &
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
    &   "WARNING::P5E - matrix-is-singular!"

#           endif

            do  ivar =  +1, nvar

            edge(ivar,ipos) =   &
        &   fdat(1,ivar,ipos-1) * 0.5d+0 + &
        &   fdat(1,ivar,ipos-0) * 0.5d+0
        
            dfdx(ivar,ipos) =   &
        &   fdat(1,ivar,ipos-0) * 0.5d+0 - &
        &   fdat(1,ivar,ipos-1) * 0.5d+0
                  
            dfdx(ivar,ipos) =   &
        &       dfdx(ivar,ipos) / xhat
        
            end do

            end if
        
            end do
        
        end if

        return
        
    end subroutine



