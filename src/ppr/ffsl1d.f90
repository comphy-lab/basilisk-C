
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
    ! FFSL1D.f90: upwind-biased flux-reconstruction scheme.
    !
    ! Darren Engwirda 
    ! 31-Mar-2019
    ! de2363 [at] columbia [dot] edu
    !
    !

    subroutine ffsl1d(npos,nvar,ndof,spac,tDEL, &
        &             mask,uvel,qbar,qedg,bclo, &
        &             bchi,work,opts)

    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell.
    ! SPAC  grid-cell spacing array. LENGTH(SPAC) == +1 if 
    !       spacing is uniform .
    ! TDEL  time-step .
    ! MASK  logical grid-cell masking array.
    ! UVEL  edge-centred velocity vectors. UVEL has SIZE = 
    !       NPOS-by-1 .
    ! QBAR  cell-centred integral moments. QBAR has SIZE =
    !       NDOF-by-NVAR-by-NPOS-1 .
    ! QEDG  edge-centred upwind flux eval. QEDG has SIZE = 
    !       NVAR-by-NPOS .
    ! BCLO  boundary condition at lower endpoint .
    ! BCHI  boundary condition at upper endpoint . 
    ! WORK  method work-space. See RCON-WORK for details .
    ! OPTS  method parameters. See RCON-OPTS for details .
    !
    
        implicit none
    
    !------------------------------------------- arguments !
        integer, intent(in)  :: npos,nvar,ndof
        class(rmap_work), intent(inout):: work
        class(rmap_opts), intent(inout):: opts
        real*8 , intent(in)  :: spac(:)
        real*8 , intent(in)  :: tDEL
        logical, intent(in)  :: mask(:)
        real*8 , intent(in)  :: qbar(:,:,:)
        real*8 , intent(in)  :: uvel(:)
        real*8 , intent(out) :: qedg(:,:)
        class(rcon_ends), intent(in) :: bclo(:)
        class(rcon_ends), intent(in) :: bchi(:)
        
    !------------------------------------------- variables !      
        integer :: head,tail,nprt
    
        head = +0 ; tail = +0 ; qedg = 0.d+0
    
        do while (.true.)

    !--------------------------------- 1. find active part !
           
            do  head = tail+1, npos-1
            if (mask(head) .eqv..true.) exit
            end do
            
            do  tail = head+1, npos-1
            if (mask(tail).neqv..true.) exit
            end do
            tail = tail - 1

            if (head.ge.npos) exit
           
    !--------------------------------- 2. rcon active part !
            
            nprt = tail - head + 1
            
            if (size(spac).ne.+1) then
            
            call rcon1d(nprt+1,nvar,ndof , &
            &    spac(    head:tail), &
            &    qbar(:,:,head:tail), &
            &    bclo,bchi,work%cell_func, &
            &    work,opts )
            
            else
            
            call rcon1d(nprt+1,nvar,ndof , &
            &    spac,qbar(:,:,head:tail), &
            &    bclo,bchi,work%cell_func, &
            &    work,opts )
            
            end if
            
    !--------------------------------- 3. int. active part !

            select case(opts%cell_meth)
                case(pcm_method)       !! 1st-order scheme
   
                if (size(spac).ne.+1) then
                
                call flux1d(nprt+1,nvar,1, &
                &    spac(  head:tail+0) , &
                &    tDEL, &
                &    uvel(  head:tail+1) , &
                &    work%cell_func, &
                &    qedg(:,head:tail+1) )
                
                else
                
                call flux1d(nprt+1,nvar,1, &
                &    spac,tDEL , &
                &    uvel(  head:tail+1) , &
                &    work%cell_func, &
                &    qedg(:,head:tail+1) )
                
                end if

                case(plm_method)       !! 2nd-order scheme
   
                if (size(spac).ne.+1) then
                
                call flux1d(nprt+1,nvar,2, &
                &    spac(  head:tail+0) , &
                &    tDEL, &
                &    uvel(  head:tail+1) , &
                &    work%cell_func, &
                &    qedg(:,head:tail+1) )
                
                else
                
                call flux1d(nprt+1,nvar,2, &
                &    spac,tDEL , &
                &    uvel(  head:tail+1) , &
                &    work%cell_func, &
                &    qedg(:,head:tail+1) )
                
                end if
       
                case(ppm_method)       !! 3rd-order scheme
   
                if (size(spac).ne.+1) then
                
                call flux1d(nprt+1,nvar,3, &
                &    spac(  head:tail+0) , &
                &    tDEL, &
                &    uvel(  head:tail+1) , &
                &    work%cell_func, &
                &    qedg(:,head:tail+1) )
                
                else
                
                call flux1d(nprt+1,nvar,3, &
                &    spac,tDEL , &
                &    uvel(  head:tail+1) , &
                &    work%cell_func, &
                &    qedg(:,head:tail+1) )
                
                end if
        
                case(pqm_method)       !! 5th-order scheme
   
                if (size(spac).ne.+1) then
                
                call flux1d(nprt+1,nvar,5, &
                &    spac(  head:tail+0) , &
                &    tDEL, &
                &    uvel(  head:tail+1) , &
                &    work%cell_func, &
                &    qedg(:,head:tail+1) )
                
                else
                
                call flux1d(nprt+1,nvar,5, &
                &    spac,tDEL , &
                &    uvel(  head:tail+1) , &
                &    work%cell_func, &
                &    qedg(:,head:tail+1) )
                
                end if
            
            end select

        end do    
  
        return
        
    end  subroutine
    
    ! FLUX1D: a degree-k, upwind-type flux reconstruction. !

    pure subroutine flux1d(npos,nvar,mdof,SPAC, &
        &                  tDEL,uvel,QHAT,qedg)
    
    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! MDOF  no. degrees-of-freedom per QHAT.
    ! SPAC  grid spacing vector. SIZE(SPAC)==+1 if uniform .
    ! TDEL  time-step .
    ! UVEL  edge-centred velocity vectors. UVEL has SIZE = 
    !       NPOS-by-1 .
    ! QHAT  cell-centred polynomial recon. QHAT has SIZE =  
    !       NDOF-by-NVAR-by-NPOS-1 .
    ! QEDG  edge-centred upwind flux eval. QEDG has SIZE = 
    !       NVAR-by-NPOS .
    !

        implicit none
        
    !------------------------------------------- arguments !
        integer, intent(in)  :: npos,nvar,mdof
        real*8 , intent(in)  :: SPAC(:)
        real*8 , intent(in)  :: tDEL
        real*8 , intent(in)  :: uvel(:)
        real*8 , intent(in)  :: QHAT(:,:,:)
        real*8 , intent(out) :: qedg(:,:)  
  
    !------------------------------------------- variables !      
        integer :: ipos,ivar
        real*8  :: uCFL,xhat,ss11,ss22,flux
        real*8  :: vv11(1:5)
        real*8  :: vv22(1:5)
        real*8  :: ivec(1:5)

    !----------- single-cell, lagrangian-type upwind rcon. !
        
        do  ipos = +2 , npos - 1
        
            if (uvel(ipos) .gt. +0.d0) then
            
    !----------- integrate profile over upwind cell IPOS-1 !
 
            if (size(SPAC).ne.+1) then
            xhat = .5d0 * SPAC(ipos-1)          
            uCFL = uvel(ipos) &
        &        * tDEL / SPAC(ipos-1)
            else
            xhat = .5d0 * SPAC(    +1)
            uCFL = uvel(ipos) &
        &        * tDEL / SPAC(    +1)
            end if
 
            ss11 = +1.d0 - 2.d0 * uCFL
            ss22 = +1.d0

            call bfun1d(-1,mdof,ss11,vv11)
            call bfun1d(-1,mdof,ss22,vv22)
            
            ivec =  vv22 - vv11
    
            do  ivar = +1, nvar

                flux =  dot_product (  &
        &           ivec(1:mdof), & 
        &       QHAT(1:mdof,ivar,ipos-1) )

                flux = flux * xhat

                qedg(ivar,ipos) = flux
                
            end do
                      
            else &
        &   if (uvel(ipos) .lt. -0.d0) then
            
    !----------- integrate profile over upwind cell IPOS+0 !
            
            if (size(SPAC).ne.+1) then
            xhat = .5d0 * SPAC(ipos-0)          
            uCFL = uvel(ipos) &
        &        * tDEL / SPAC(ipos-0)
            else
            xhat = .5d0 * SPAC(    +1)
            uCFL = uvel(ipos) &
        &        * tDEL / SPAC(    +1)
            end if

            ss11 = -1.d0 - 2.d0 * uCFL
            ss22 = -1.d0
        
            call bfun1d(-1,mdof,ss11,vv11)
            call bfun1d(-1,mdof,ss22,vv22)
            
            ivec =  vv22 - vv11
    
            do  ivar = +1, nvar

                flux =  dot_product (  &
        &           ivec(1:mdof), & 
        &       QHAT(1:mdof,ivar,ipos-0) )

                flux = flux * xhat

                qedg(ivar,ipos) = flux
                
            end do
            
            end if
        
        end do
               
        return

    end  subroutine



