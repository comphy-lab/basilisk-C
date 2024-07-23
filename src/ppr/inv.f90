
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
    ! INV.f90: block-wise solution of small linear systems.
    !
    ! Darren Engwirda 
    ! 25-Mar-2019
    ! de2363 [at] columbia [dot] edu
    !
    !

    pure subroutine inv_2x2(amat,adim,ainv,vdim, &
        &                   adet)

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: adim
        real*8 , intent( in) :: amat(adim,*)
        integer, intent( in) :: vdim
        real*8 , intent(out) :: ainv(vdim,*)
        real*8 , intent(out) :: adet
        
    !------------------------------------------- form A^-1 !

        adet   = amat(1,1) * amat(2,2) &
               - amat(1,2) * amat(2,1)

        ainv(1,1) =          amat(2,2)
        ainv(1,2) =        - amat(1,2)
        ainv(2,1) =        - amat(2,1)
        ainv(2,2) =          amat(1,1)

        return

    end  subroutine

    pure subroutine inv_3x3(amat,adim,ainv,vdim, &
        &                   adet)

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: adim
        real*8 , intent( in) :: amat(adim,*)
        integer, intent( in) :: vdim
        real*8 , intent(out) :: ainv(vdim,*)
        real*8 , intent(out) :: adet
        
    !------------------------------------------- variables !
        real*8 :: &
        aa2233,aa2332,aa2133,aa2331,aa2132,&
        aa2231,aa1233,aa1332,aa1223,aa1322,&
        aa1133,aa1331,aa1123,aa1321,aa1132,&
        aa1231,aa1122,aa1221

    !------------------------------------------- form A^-1 !

        aa2233 = amat(2,2) * amat(3,3)
        aa2332 = amat(2,3) * amat(3,2)
        aa2133 = amat(2,1) * amat(3,3)
        aa2331 = amat(2,3) * amat(3,1)
        aa2132 = amat(2,1) * amat(3,2)
        aa2231 = amat(2,2) * amat(3,1)

        adet =  &
        amat(1,1) *  (aa2233 - aa2332) - &
	    amat(1,2) *  (aa2133 - aa2331) + &
	    amat(1,3) *  (aa2132 - aa2231)

        aa1233 = amat(1,2) * amat(3,3)
        aa1332 = amat(1,3) * amat(3,2)
        aa1223 = amat(1,2) * amat(2,3)
        aa1322 = amat(1,3) * amat(2,2)
        aa1133 = amat(1,1) * amat(3,3)
        aa1331 = amat(1,3) * amat(3,1)
        aa1123 = amat(1,1) * amat(2,3)
        aa1321 = amat(1,3) * amat(2,1)
        aa1132 = amat(1,1) * amat(3,2)
        aa1231 = amat(1,2) * amat(3,1)
        aa1122 = amat(1,1) * amat(2,2)
        aa1221 = amat(1,2) * amat(2,1)

        ainv(1,1) =  (aa2233 - aa2332)
        ainv(1,2) = -(aa1233 - aa1332)
        ainv(1,3) =  (aa1223 - aa1322)
        
        ainv(2,1) = -(aa2133 - aa2331)
        ainv(2,2) =  (aa1133 - aa1331)
        ainv(2,3) = -(aa1123 - aa1321)
        
        ainv(3,1) =  (aa2132 - aa2231)
        ainv(3,2) = -(aa1132 - aa1231)
        ainv(3,3) =  (aa1122 - aa1221)

        return

    end  subroutine

    pure subroutine mul_2x2(amat,adim,bmat,bdim, &
        &                   scal,cmat,cdim)

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)    :: adim
        real*8 , intent(in)    :: amat(adim,*)
        integer, intent(in)    :: bdim
        real*8 , intent(in)    :: bmat(bdim,*)
        real*8 , intent(in)    :: scal
        integer, intent(in)    :: cdim
        real*8 , intent(inout) :: cmat(cdim,*)

    !-------------------------------- C = C + scal * A * B !

        if (scal .eq. +1.d0) then

        cmat(1,1) = cmat(1,1) & 
             + ( amat(1,1) * bmat(1,1) &
               + amat(1,2) * bmat(2,1) )
        cmat(2,1) = cmat(2,1) & 
             + ( amat(2,1) * bmat(1,1) &
               + amat(2,2) * bmat(2,1) )

        cmat(1,2) = cmat(1,2) & 
             + ( amat(1,1) * bmat(1,2) &
               + amat(1,2) * bmat(2,2) )
        cmat(2,2) = cmat(2,2) & 
             + ( amat(2,1) * bmat(1,2) &
               + amat(2,2) * bmat(2,2) )

        else &
        if (scal .eq. -1.d0) then

        cmat(1,1) = cmat(1,1) & 
             - ( amat(1,1) * bmat(1,1) &
               + amat(1,2) * bmat(2,1) )
        cmat(2,1) = cmat(2,1) & 
             - ( amat(2,1) * bmat(1,1) &
               + amat(2,2) * bmat(2,1) )

        cmat(1,2) = cmat(1,2) & 
             - ( amat(1,1) * bmat(1,2) &
               + amat(1,2) * bmat(2,2) )
        cmat(2,2) = cmat(2,2) & 
             - ( amat(2,1) * bmat(1,2) &
               + amat(2,2) * bmat(2,2) )

        else 

        cmat(1,1) = cmat(1,1) + & 
        scal * ( amat(1,1) * bmat(1,1) &
               + amat(1,2) * bmat(2,1) )
        cmat(2,1) = cmat(2,1) + & 
        scal * ( amat(2,1) * bmat(1,1) &
               + amat(2,2) * bmat(2,1) )

        cmat(1,2) = cmat(1,2) + & 
        scal * ( amat(1,1) * bmat(1,2) &
               + amat(1,2) * bmat(2,2) )
        cmat(2,2) = cmat(2,2) + & 
        scal * ( amat(2,1) * bmat(1,2) &
               + amat(2,2) * bmat(2,2) )

        end if

        return

    end  subroutine

    pure subroutine mul_3x3(amat,adim,bmat,bdim, &
        &                   scal,cmat,cdim)

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)    :: adim
        real*8 , intent(in)    :: amat(adim,*)
        integer, intent(in)    :: bdim
        real*8 , intent(in)    :: bmat(bdim,*)
        real*8 , intent(in)    :: scal
        integer, intent(in)    :: cdim
        real*8 , intent(inout) :: cmat(cdim,*)

    !-------------------------------- C = C + scal * A * B !

        if (scal .eq. +1.d0) then

        cmat(1,1) = cmat(1,1) &
             + ( amat(1,1) * bmat(1,1) &
               + amat(1,2) * bmat(2,1) &
               + amat(1,3) * bmat(3,1) )
        cmat(2,1) = cmat(2,1) &
             + ( amat(2,1) * bmat(1,1) &
               + amat(2,2) * bmat(2,1) &
               + amat(2,3) * bmat(3,1) )
        cmat(3,1) = cmat(3,1) &
             + ( amat(3,1) * bmat(1,1) &
               + amat(3,2) * bmat(2,1) &
               + amat(3,3) * bmat(3,1) )

        cmat(1,2) = cmat(1,2) &
             + ( amat(1,1) * bmat(1,2) &
               + amat(1,2) * bmat(2,2) &
               + amat(1,3) * bmat(3,2) )
        cmat(2,2) = cmat(2,2) &
             + ( amat(2,1) * bmat(1,2) &
               + amat(2,2) * bmat(2,2) &
               + amat(2,3) * bmat(3,2) )
        cmat(3,2) = cmat(3,2) &
             + ( amat(3,1) * bmat(1,2) &
               + amat(3,2) * bmat(2,2) &
               + amat(3,3) * bmat(3,2) )

        cmat(1,3) = cmat(1,3) &
             + ( amat(1,1) * bmat(1,3) &
               + amat(1,2) * bmat(2,3) &
               + amat(1,3) * bmat(3,3) )
        cmat(2,3) = cmat(2,3) &
             + ( amat(2,1) * bmat(1,3) &
               + amat(2,2) * bmat(2,3) &
               + amat(2,3) * bmat(3,3) )
        cmat(3,3) = cmat(3,3) &
             + ( amat(3,1) * bmat(1,3) &
               + amat(3,2) * bmat(2,3) &
               + amat(3,3) * bmat(3,3) )

        else &
        if (scal .eq. -1.d0) then

        cmat(1,1) = cmat(1,1) &
             - ( amat(1,1) * bmat(1,1) &
               + amat(1,2) * bmat(2,1) &
               + amat(1,3) * bmat(3,1) )
        cmat(2,1) = cmat(2,1) &
             - ( amat(2,1) * bmat(1,1) &
               + amat(2,2) * bmat(2,1) &
               + amat(2,3) * bmat(3,1) )
        cmat(3,1) = cmat(3,1) &
             - ( amat(3,1) * bmat(1,1) &
               + amat(3,2) * bmat(2,1) &
               + amat(3,3) * bmat(3,1) )

        cmat(1,2) = cmat(1,2) &
             - ( amat(1,1) * bmat(1,2) &
               + amat(1,2) * bmat(2,2) &
               + amat(1,3) * bmat(3,2) )
        cmat(2,2) = cmat(2,2) &
             - ( amat(2,1) * bmat(1,2) &
               + amat(2,2) * bmat(2,2) &
               + amat(2,3) * bmat(3,2) )
        cmat(3,2) = cmat(3,2) &
             - ( amat(3,1) * bmat(1,2) &
               + amat(3,2) * bmat(2,2) &
               + amat(3,3) * bmat(3,2) )

        cmat(1,3) = cmat(1,3) &
             - ( amat(1,1) * bmat(1,3) &
               + amat(1,2) * bmat(2,3) &
               + amat(1,3) * bmat(3,3) )
        cmat(2,3) = cmat(2,3) &
             - ( amat(2,1) * bmat(1,3) &
               + amat(2,2) * bmat(2,3) &
               + amat(2,3) * bmat(3,3) )
        cmat(3,3) = cmat(3,3) &
             - ( amat(3,1) * bmat(1,3) &
               + amat(3,2) * bmat(2,3) &
               + amat(3,3) * bmat(3,3) )

        else 

        cmat(1,1) = cmat(1,1) + & 
        scal * ( amat(1,1) * bmat(1,1) &
               + amat(1,2) * bmat(2,1) &
               + amat(1,3) * bmat(3,1) )
        cmat(2,1) = cmat(2,1) + & 
        scal * ( amat(2,1) * bmat(1,1) &
               + amat(2,2) * bmat(2,1) &
               + amat(2,3) * bmat(3,1) )
        cmat(3,1) = cmat(3,1) + & 
        scal * ( amat(3,1) * bmat(1,1) &
               + amat(3,2) * bmat(2,1) &
               + amat(3,3) * bmat(3,1) )

        cmat(1,2) = cmat(1,2) + & 
        scal * ( amat(1,1) * bmat(1,2) &
               + amat(1,2) * bmat(2,2) &
               + amat(1,3) * bmat(3,2) )
        cmat(2,2) = cmat(2,2) + & 
        scal * ( amat(2,1) * bmat(1,2) &
               + amat(2,2) * bmat(2,2) &
               + amat(2,3) * bmat(3,2) )
        cmat(3,2) = cmat(3,2) + & 
        scal * ( amat(3,1) * bmat(1,2) &
               + amat(3,2) * bmat(2,2) &
               + amat(3,3) * bmat(3,2) )

        cmat(1,3) = cmat(1,3) + & 
        scal * ( amat(1,1) * bmat(1,3) &
               + amat(1,2) * bmat(2,3) &
               + amat(1,3) * bmat(3,3) )
        cmat(2,3) = cmat(2,3) + & 
        scal * ( amat(2,1) * bmat(1,3) &
               + amat(2,2) * bmat(2,3) &
               + amat(2,3) * bmat(3,3) )
        cmat(3,3) = cmat(3,3) + & 
        scal * ( amat(3,1) * bmat(1,3) &
               + amat(3,2) * bmat(2,3) &
               + amat(3,3) * bmat(3,3) )

        end if

        return

    end  subroutine

    pure subroutine slv_2x2(amat,adim,vrhs,vdim, & 
        &                   nrhs,fEPS,okay)

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)    :: adim
        real*8 , intent(in)    :: amat(adim,*)
        integer, intent(in)    :: vdim
        real*8 , intent(inout) :: vrhs(vdim,*)
        integer, intent(in)    :: nrhs
        real*8 , intent(in)    :: fEPS
        logical, intent(inout) :: okay

    !------------------------------------------- variables !
        real*8                 :: ainv(2,2)
        real*8                 :: adet
        real*8                 :: vtmp(  2)
        integer                :: irhs

        integer, parameter     :: LDIM = 2

    !---------------------------------------- calc. inv(A) !

        call inv_2x2(amat,adim,ainv,LDIM,&
                     adet)

        okay = (abs(adet) .gt. fEPS)
        
        if (okay.eqv..false.) return
    
    !---------------------------------------- v = A^-1 * v !
    
        do irhs = 1, nrhs

        vtmp(1) =  &
          + ( & 
        ainv(1, 1) *   vrhs(1,irhs) &
      + ainv(1, 2) *   vrhs(2,irhs) &
            ) / adet

        vtmp(2) =  &
          + ( & 
        ainv(2, 1) *   vrhs(1,irhs) &
      + ainv(2, 2) *   vrhs(2,irhs) &
            ) / adet

        vrhs(1,irhs) = vtmp(1)
        vrhs(2,irhs) = vtmp(2)

        end do

        return

    end  subroutine

    pure subroutine slv_3x3(amat,adim,vrhs,vdim, & 
        &                   nrhs,fEPS,okay)

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)    :: adim
        real*8 , intent(in)    :: amat(adim,*)
        integer, intent(in)    :: vdim
        real*8 , intent(inout) :: vrhs(vdim,*)
        integer, intent(in)    :: nrhs
        real*8 , intent(in)    :: fEPS
        logical, intent(inout) :: okay

    !------------------------------------------- variables !
        real*8                 :: ainv(3,3)
        real*8                 :: adet
        real*8                 :: vtmp(  3)
        integer                :: irhs

        integer, parameter     :: LDIM = 3

    !---------------------------------------- calc. inv(A) !

        call inv_3x3(amat,adim,ainv,LDIM,&
                     adet)

        okay = (abs(adet) .gt. fEPS)
        
        if (okay.eqv..false.) return
    
    !---------------------------------------- v = A^-1 * v !
    
        do irhs = 1, nrhs

        vtmp(1) =  &
          + ( & 
        ainv(1, 1) *   vrhs(1,irhs) &
      + ainv(1, 2) *   vrhs(2,irhs) &
      + ainv(1, 3) *   vrhs(3,irhs) &
            ) / adet

        vtmp(2) =  &
          + ( & 
        ainv(2, 1) *   vrhs(1,irhs) &
      + ainv(2, 2) *   vrhs(2,irhs) &
      + ainv(2, 3) *   vrhs(3,irhs) &
            ) / adet

        vtmp(3) =  &
          + ( & 
        ainv(3, 1) *   vrhs(1,irhs) &
      + ainv(3, 2) *   vrhs(2,irhs) &
      + ainv(3, 3) *   vrhs(3,irhs) &
            ) / adet

        vrhs(1,irhs) = vtmp(1)
        vrhs(2,irhs) = vtmp(2)
        vrhs(3,irhs) = vtmp(3)

        end do

        return

    end  subroutine

    pure subroutine slv_4x4(amat,adim,vrhs,vdim, & 
        &                   nrhs,fEPS,okay)

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)    :: adim
        real*8 , intent(in)    :: amat(adim,*)
        integer, intent(in)    :: vdim
        real*8 , intent(inout) :: vrhs(vdim,*)
        integer, intent(in)    :: nrhs
        real*8 , intent(in)    :: fEPS
        logical, intent(inout) :: okay

    !------------------------------------------- variables !
        real*8                 :: ainv(2,2)
        real*8                 :: lmat(2,2)
        real*8                 :: umat(2,2)
        real*8                 :: smat(2,2)
        real*8                 :: sinv(2,2)
        real*8                 :: adet,sdet        
        real*8                 :: vtmp(  2)
        integer                :: irhs

        integer, parameter     :: LDIM = 2

    !---------------------- form a block LDU factorisation !

        call inv_2x2(amat(1,1),adim,ainv,LDIM, &
                     adet)

        okay = (abs(adet) .gt. fEPS)
        
        if (okay.eqv..false.) return

    !---------------------------------------- L = C * A^-1 !

        lmat(1,1) = +0.d0
        lmat(1,2) = +0.d0
        lmat(2,1) = +0.d0
        lmat(2,2) = +0.d0

        call mul_2x2(amat(3,1),adim,ainv,LDIM, &
                    +1.d0,lmat,LDIM)
        
    !---------------------------------------- U = A^-1 * B !       

        umat(1,1) = +0.d0
        umat(1,2) = +0.d0
        umat(2,1) = +0.d0
        umat(2,2) = +0.d0
 
        call mul_2x2(ainv,LDIM,amat(1,3),adim, &
                    +1.d0,umat,LDIM)
        
    !-------------------------------- S = D - C * A^-1 * B !

        smat(1,1) = amat(3,3)
        smat(1,2) = amat(3,4)
        smat(2,1) = amat(4,3)
        smat(2,2) = amat(4,4)

        call mul_2x2(lmat,LDIM,amat(1,3),adim, &
                    -1.d0/adet,smat,LDIM)

        call inv_2x2(smat,LDIM,sinv,LDIM,sdet)

        okay = (abs(adet) .gt. fEPS)
        
        if (okay.eqv..false.) return

    !-------------------------------- back-solve LDU = rhs !

        do irhs = 1, nrhs

    !---------------------------------------- solve L part !

        vrhs(3,irhs) = vrhs(3,irhs) &
          - ( &
        lmat(1, 1) *   vrhs(1,irhs) &
      + lmat(1, 2) *   vrhs(2,irhs) &
            ) / adet

        vrhs(4,irhs) = vrhs(4,irhs) &
          - ( &
        lmat(2, 1) *   vrhs(1,irhs) &
      + lmat(2, 2) *   vrhs(2,irhs) &
            ) / adet
        
    !---------------------------------------- solve D part !

        vtmp(1) =  &
          + ( & 
        ainv(1, 1) *   vrhs(1,irhs) &
      + ainv(1, 2) *   vrhs(2,irhs) &
            ) / adet

        vtmp(2) =  &
          + ( & 
        ainv(2, 1) *   vrhs(1,irhs) &
      + ainv(2, 2) *   vrhs(2,irhs) &
            ) / adet

        vrhs(1,irhs) = vtmp(1)
        vrhs(2,irhs) = vtmp(2)

        vtmp(1) =  &
          + ( &
        sinv(1, 1) *   vrhs(3,irhs) &
      + sinv(1, 2) *   vrhs(4,irhs) & 
            ) / sdet

        vtmp(2) =  &
          + ( &
        sinv(2, 1) *   vrhs(3,irhs) &
      + sinv(2, 2) *   vrhs(4,irhs) & 
            ) / sdet

        vrhs(3,irhs) = vtmp(1)
        vrhs(4,irhs) = vtmp(2)

    !---------------------------------------- solve U part !

        vrhs(1,irhs) = vrhs(1,irhs) &
          - ( &
        umat(1, 1) *   vrhs(3,irhs) &
      + umat(1, 2) *   vrhs(4,irhs) & 
            ) / adet

        vrhs(2,irhs) = vrhs(2,irhs) &
          - ( &
        umat(2, 1) *   vrhs(3,irhs) &
      + umat(2, 2) *   vrhs(4,irhs) & 
            ) / adet

        end do

        return

    end  subroutine

    pure subroutine slv_6x6(amat,adim,vrhs,vdim, & 
        &                   nrhs,fEPS,okay)

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)    :: adim
        real*8 , intent(in)    :: amat(adim,*)
        integer, intent(in)    :: vdim
        real*8 , intent(inout) :: vrhs(vdim,*)
        integer, intent(in)    :: nrhs
        real*8 , intent(in)    :: fEPS
        logical, intent(inout) :: okay

    !------------------------------------------- variables !
        real*8                 :: ainv(3,3)
        real*8                 :: lmat(3,3)
        real*8                 :: umat(3,3)
        real*8                 :: smat(3,3)
        real*8                 :: sinv(3,3)
        real*8                 :: adet,sdet        
        real*8                 :: vtmp(  3)
        integer                :: irhs

        integer, parameter     :: LDIM = 3

    !---------------------- form a block LDU factorisation !

        call inv_3x3(amat(1,1),adim,ainv,LDIM, &
                     adet)

        okay = (abs(adet) .gt. fEPS)
        
        if (okay.eqv..false.) return

    !---------------------------------------- L = C * A^-1 !

        lmat(1,1) = +0.d0
        lmat(1,2) = +0.d0
        lmat(1,3) = +0.d0
        lmat(2,1) = +0.d0
        lmat(2,2) = +0.d0
        lmat(2,3) = +0.d0
        lmat(3,1) = +0.d0
        lmat(3,2) = +0.d0
        lmat(3,3) = +0.d0

        call mul_3x3(amat(4,1),adim,ainv,LDIM, &
                    +1.d0,lmat,LDIM)
        
    !---------------------------------------- U = A^-1 * B !       

        umat(1,1) = +0.d0
        umat(1,2) = +0.d0
        umat(1,3) = +0.d0
        umat(2,1) = +0.d0
        umat(2,2) = +0.d0
        umat(2,3) = +0.d0
        umat(3,1) = +0.d0
        umat(3,2) = +0.d0
        umat(3,3) = +0.d0
 
        call mul_3x3(ainv,LDIM,amat(1,4),adim, &
                    +1.d0,umat,LDIM)
        
    !-------------------------------- S = D - C * A^-1 * B !

        smat(1,1) = amat(4,4)
        smat(1,2) = amat(4,5)
        smat(1,3) = amat(4,6)
        smat(2,1) = amat(5,4)
        smat(2,2) = amat(5,5)
        smat(2,3) = amat(5,6)
        smat(3,1) = amat(6,4)
        smat(3,2) = amat(6,5)
        smat(3,3) = amat(6,6)

        call mul_3x3(lmat,LDIM,amat(1,4),adim, &
                    -1.d0/adet,smat,LDIM)

        call inv_3x3(smat,LDIM,sinv,LDIM,sdet)

        okay = (abs(adet) .gt. fEPS)
        
        if (okay.eqv..false.) return

    !-------------------------------- back-solve LDU = rhs !

        do irhs = 1, nrhs

    !---------------------------------------- solve L part !

        vrhs(4,irhs) = vrhs(4,irhs) &
          - ( &
        lmat(1, 1) *   vrhs(1,irhs) &
      + lmat(1, 2) *   vrhs(2,irhs) &
      + lmat(1, 3) *   vrhs(3,irhs) &
            ) / adet

        vrhs(5,irhs) = vrhs(5,irhs) &
          - ( &
        lmat(2, 1) *   vrhs(1,irhs) &
      + lmat(2, 2) *   vrhs(2,irhs) &
      + lmat(2, 3) *   vrhs(3,irhs) &
            ) / adet

        vrhs(6,irhs) = vrhs(6,irhs) &
          - ( &
        lmat(3, 1) *   vrhs(1,irhs) &
      + lmat(3, 2) *   vrhs(2,irhs) &
      + lmat(3, 3) *   vrhs(3,irhs) &
            ) / adet
        
    !---------------------------------------- solve D part !

        vtmp(1) =  &
          + ( & 
        ainv(1, 1) *   vrhs(1,irhs) &
      + ainv(1, 2) *   vrhs(2,irhs) &
      + ainv(1, 3) *   vrhs(3,irhs) &
            ) / adet

        vtmp(2) =  &
          + ( & 
        ainv(2, 1) *   vrhs(1,irhs) &
      + ainv(2, 2) *   vrhs(2,irhs) &
      + ainv(2, 3) *   vrhs(3,irhs) &
            ) / adet

        vtmp(3) =  &
          + ( & 
        ainv(3, 1) *   vrhs(1,irhs) &
      + ainv(3, 2) *   vrhs(2,irhs) &
      + ainv(3, 3) *   vrhs(3,irhs) &
            ) / adet

        vrhs(1,irhs) = vtmp(1)
        vrhs(2,irhs) = vtmp(2)
        vrhs(3,irhs) = vtmp(3)

        vtmp(1) =  &
          + ( &
        sinv(1, 1) *   vrhs(4,irhs) &
      + sinv(1, 2) *   vrhs(5,irhs) &
      + sinv(1, 3) *   vrhs(6,irhs) & 
            ) / sdet

        vtmp(2) =  &
          + ( &
        sinv(2, 1) *   vrhs(4,irhs) &
      + sinv(2, 2) *   vrhs(5,irhs) &
      + sinv(2, 3) *   vrhs(6,irhs) & 
            ) / sdet

        vtmp(3) =  &
          + ( &
        sinv(3, 1) *   vrhs(4,irhs) &
      + sinv(3, 2) *   vrhs(5,irhs) &
      + sinv(3, 3) *   vrhs(6,irhs) & 
            ) / sdet

        vrhs(4,irhs) = vtmp(1)
        vrhs(5,irhs) = vtmp(2)
        vrhs(6,irhs) = vtmp(3)

    !---------------------------------------- solve U part !

        vrhs(1,irhs) = vrhs(1,irhs) &
          - ( &
        umat(1, 1) *   vrhs(4,irhs) &
      + umat(1, 2) *   vrhs(5,irhs) &
      + umat(1, 3) *   vrhs(6,irhs) & 
            ) / adet

        vrhs(2,irhs) = vrhs(2,irhs) &
          - ( &
        umat(2, 1) *   vrhs(4,irhs) &
      + umat(2, 2) *   vrhs(5,irhs) &
      + umat(2, 3) *   vrhs(6,irhs) & 
            ) / adet

        vrhs(3,irhs) = vrhs(3,irhs) &
          - ( &
        umat(3, 1) *   vrhs(4,irhs) &
      + umat(3, 2) *   vrhs(5,irhs) &
      + umat(3, 3) *   vrhs(6,irhs) & 
            ) / adet

        end do

        return

    end  subroutine



