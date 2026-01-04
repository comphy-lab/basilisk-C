
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
    ! ROOT1D.f90: find the "roots" of degree-k polynomials.
    !
    ! Darren Engwirda 
    ! 25-Mar-2019
    ! de2363 [at] columbia [dot] edu
    !
    !

    pure subroutine roots_2(aa,bb,cc,xx,haveroot)

    !
    ! solve:: aa * xx**2 + bb * xx**1 + cc = +0.0 .
    !

        implicit none

    !------------------------------------------- arguments !
        real*8 , intent( in) :: aa,bb,cc
        real*8 , intent(out) :: xx(1:2)
        logical, intent(out) :: haveroot

    !------------------------------------------- variables !
        real*8 :: sq,ia,a0,b0,c0,x0

        real*8, parameter :: rt = +1.d-14

        a0 = abs(aa)
        b0 = abs(bb)
        c0 = abs(cc)

        sq = bb * bb - 4.0d+0 * aa * cc

        if (sq .ge. 0.0d+0) then

            sq = sqrt (sq)

            xx(1) =  - bb + sq
            xx(2) =  - bb - sq

            x0 = max(abs(xx(1)), &
        &            abs(xx(2)))

            if (a0 .gt. (rt*x0)) then

    !-------------------------------------- degree-2 roots !
    
            haveroot =  .true.

            ia = 0.5d+0   / aa

            xx(1) = xx(1) * ia
            xx(2) = xx(2) * ia
            
            else &
        &   if (b0 .gt. (rt*c0)) then

    !-------------------------------------- degree-1 roots !

            haveroot =  .true.
            
            xx(1) =  - cc / bb
            xx(2) =  - cc / bb
            
            else
            
            haveroot = .false.
            
            end if

        else

            haveroot = .false.

        end if

        return

    end subroutine



