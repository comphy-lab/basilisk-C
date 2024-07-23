
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
    ! BFUN1D.f90: poly. basis-functions for reconstruction.
    !
    ! Darren Engwirda 
    ! 07-Sep-2016
    ! de2363 [at] columbia [dot] edu
    !
    !

    pure subroutine bfun1d(isel,ndof,sval,bfun)

    !
    ! ISEL  basis-function "order", -1 => integral-basis , 
    !       +0 => function-basis, +1 => 1st deriv.-basis ,
    !       +2 => 2nd deriv.-basis.
    ! NDOF  no. degrees-of-freedom in basis.
    ! SVAL  local coord. at which to evaluate basis-func.,
    !       such that -1.0 <= SVAL <= +1.0 .
    ! BFUN  basis-vector evaluated at SVAL .
    !
    
        implicit none
        
    !------------------------------------------- arguments !
        integer, intent( in) :: isel,ndof
        real*8 , intent( in) :: sval
        real*8 , intent(out) :: bfun(:)
        
        select case (isel)
        case (-1)
    !------------------------------------ -1th-order basis !
            select case (ndof)
            case (+1)
                bfun(1) = sval**1 / 1.d0
                
            case (+2)
                bfun(1) = sval**1 / 1.d0
                bfun(2) = sval**2 / 2.d0
                
            case (+3)
                bfun(1) = sval**1 / 1.d0
                bfun(2) = sval**2 / 2.d0
                bfun(3) = sval**3 / 3.d0
                
            case (+4)
                bfun(1) = sval**1 / 1.d0
                bfun(2) = sval**2 / 2.d0
                bfun(3) = sval**3 / 3.d0
                bfun(4) = sval**4 / 4.d0
                
            case (+5)
                bfun(1) = sval**1 / 1.d0
                bfun(2) = sval**2 / 2.d0
                bfun(3) = sval**3 / 3.d0
                bfun(4) = sval**4 / 4.d0
                bfun(5) = sval**5 / 5.d0
                
            case (+6)
                bfun(1) = sval**1 / 1.d0
                bfun(2) = sval**2 / 2.d0
                bfun(3) = sval**3 / 3.d0
                bfun(4) = sval**4 / 4.d0
                bfun(5) = sval**5 / 5.d0
                bfun(6) = sval**6 / 6.d0
                
            case (+7)
                bfun(1) = sval**1 / 1.d0
                bfun(2) = sval**2 / 2.d0
                bfun(3) = sval**3 / 3.d0
                bfun(4) = sval**4 / 4.d0
                bfun(5) = sval**5 / 5.d0
                bfun(6) = sval**6 / 6.d0
                bfun(7) = sval**7 / 7.d0
            
            end select

        case (+0)
    !------------------------------------ +0th-order basis !
            select case (ndof)
            case (+1)
                bfun(1) =           1.d0
                
            case (+2)
                bfun(1) =           1.d0
                bfun(2) = sval**1 * 1.d0
                
            case (+3)
                bfun(1) =           1.d0
                bfun(2) = sval**1 * 1.d0
                bfun(3) = sval**2 * 1.d0
                
            case (+4)
                bfun(1) =           1.d0
                bfun(2) = sval**1 * 1.d0
                bfun(3) = sval**2 * 1.d0
                bfun(4) = sval**3 * 1.d0
                
            case (+5)
                bfun(1) =           1.d0
                bfun(2) = sval**1 * 1.d0
                bfun(3) = sval**2 * 1.d0
                bfun(4) = sval**3 * 1.d0
                bfun(5) = sval**4 * 1.d0
                
            case (+6)
                bfun(1) =           1.d0
                bfun(2) = sval**1 * 1.d0
                bfun(3) = sval**2 * 1.d0
                bfun(4) = sval**3 * 1.d0
                bfun(5) = sval**4 * 1.d0
                bfun(6) = sval**5 * 1.d0
                
            case (+7)
                bfun(1) =           1.d0
                bfun(2) = sval**1 * 1.d0
                bfun(3) = sval**2 * 1.d0
                bfun(4) = sval**3 * 1.d0
                bfun(5) = sval**4 * 1.d0
                bfun(6) = sval**5 * 1.d0
                bfun(7) = sval**6 * 1.d0
            
            end select

        case (+1)
    !------------------------------------ +1st-order basis !
            select case (ndof)
            case (+1)
                bfun(1) =           0.d0
                
            case (+2)
                bfun(1) =           0.d0
                bfun(2) =           1.d0
                
            case (+3)
                bfun(1) =           0.d0
                bfun(2) =           1.d0
                bfun(3) = sval**1 * 2.d0
                
            case (+4)
                bfun(1) =           0.d0
                bfun(2) =           1.d0
                bfun(3) = sval**1 * 2.d0
                bfun(4) = sval**2 * 3.d0
                
            case (+5)
                bfun(1) =           0.d0
                bfun(2) =           1.d0
                bfun(3) = sval**1 * 2.d0
                bfun(4) = sval**2 * 3.d0
                bfun(5) = sval**3 * 4.d0
                
            case (+6)
                bfun(1) =           0.d0
                bfun(2) =           1.d0
                bfun(3) = sval**1 * 2.d0
                bfun(4) = sval**2 * 3.d0
                bfun(5) = sval**3 * 4.d0
                bfun(6) = sval**4 * 5.d0
                
            case (+7)
                bfun(1) =           0.d0
                bfun(2) =           1.d0
                bfun(3) = sval**1 * 2.d0
                bfun(4) = sval**2 * 3.d0
                bfun(5) = sval**3 * 4.d0
                bfun(6) = sval**4 * 5.d0
                bfun(7) = sval**5 * 6.d0
            
            end select

        case (+2)
    !------------------------------------ +2nd-order basis !
            select case (ndof)
            case (+1)
                bfun(1) =           0.d0
                
            case (+2)
                bfun(1) =           0.d0
                bfun(2) =           0.d0
                
            case (+3)
                bfun(1) =           0.d0
                bfun(2) =           0.d0
                bfun(3) =           2.d0
                
            case (+4)
                bfun(1) =           0.d0
                bfun(2) =           0.d0
                bfun(3) =           2.d0
                bfun(4) = sval**1 * 6.d0
                
            case (+5)
                bfun(1) =           0.d0
                bfun(2) =           0.d0
                bfun(3) =           2.d0
                bfun(4) = sval**1 * 6.d0
                bfun(5) = sval**2 *12.d0
                
            case (+6)
                bfun(1) =           0.d0
                bfun(2) =           0.d0
                bfun(3) =           2.d0
                bfun(4) = sval**1 * 6.d0
                bfun(5) = sval**2 *12.d0
                bfun(6) = sval**3 *20.d0
                
            case (+7)
                bfun(1) =           0.d0
                bfun(2) =           0.d0
                bfun(3) =           2.d0
                bfun(4) = sval**1 * 6.d0
                bfun(5) = sval**2 *12.d0
                bfun(6) = sval**3 *20.d0
                bfun(7) = sval**4 *30.d0
            
            end select

        end select
    
    end subroutine
    
    
    
