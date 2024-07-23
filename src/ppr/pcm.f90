
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
    ! PCM.f90: 1d piecewise constant reconstruction .
    !
    ! Darren Engwirda 
    ! 08-Sep-2016
    ! de2363 [at] columbia [dot] edu
    !
    !

    pure subroutine pcm(npos,nvar,ndof,fdat, &
        &               fhat)

    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! FHAT  grid-cell re-con. array. FHAT is an array with
    !       SIZE = MDOF-by-NVAR-by-NPOS-1 .
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        real*8 , intent(out) :: fhat(:,:,:)
        real*8 , intent( in) :: fdat(:,:,:)

    !------------------------------------------- variables !
        integer:: ipos,ivar,idof

        do  ipos = +1, npos - 1
        do  ivar = +1, nvar + 0
        do  idof = +1, ndof + 0

            fhat(idof,ivar,ipos) = fdat(idof,ivar,ipos)

        end do
        end do
        end do
        
        return

    end subroutine
    
    
    
