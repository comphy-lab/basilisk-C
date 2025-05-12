
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
    ! RMAP1D.f90: high-order integral re-mapping operators.
    !
    ! Darren Engwirda 
    ! 31-Mar-2019
    ! ​de2363 [at] columbia [dot] edu
    !
    !

    subroutine rmap1d(npos,nnew,nvar,ndof,xpos, &
        &             xnew,fdat,fnew,bclo,bcup, &
        &             work,opts,tCPU)

    !
    ! NPOS  no. edges in old grid.
    ! NNEW  no. edges in new grid.
    ! NVAR  no. discrete variables to remap.
    ! NDOF  no. degrees-of-freedom per cell.
    ! XPOS  old grid edge positions. XPOS is a length NPOS
    !       array .
    ! XNEW  new grid edge positions. XNEW is a length NNEW
    !       array .
    ! FDAT  grid-cell moments on old grid. FNEW has SIZE = 
    !       NDOF-by-NVAR-by-NNEW-1 .
    ! FNEW  grid-cell moments on new grid. FNEW has SIZE = 
    !       NDOF-by-NVAR-by-NNEW-1 .
    ! BCLO  boundary condition at lower endpoint .
    ! BCHI  boundary condition at upper endpoint . 
    ! WORK  method work-space. See RCON-WORK for details .
    ! OPTS  method parameters. See RCON-OPTS for details .
    ! TCPU  method tcpu-timer.
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nnew
        integer, intent( in) :: nvar,ndof
        class(rmap_work), intent(inout):: work
        class(rmap_opts), intent(inout):: opts
        real*8 , intent( in) :: xpos(:)
        real*8 , intent( in) :: xnew(:)
        real*8 , intent( in) :: fdat(:,:,:)
        real*8 , intent(out) :: fnew(:,:,:)
        type (rcon_ends), intent(in) :: bclo(:)
        type (rcon_ends), intent(in) :: bcup(:)
        type (rmap_tics), &
        &   intent(inout) , optional :: tCPU

        real*8 , parameter :: RTOL = +1.d-14

    !------------------------------------------- variables !
        integer :: ipos
        real*8  :: diff,spac,same,xtol
        real*8  :: delx(1)
        logical :: uniform
        
#       ifdef __PPR_TIMER__
        integer(kind=8) :: ttic,ttoc,rate
#       endif

        if (ndof.lt.1) return
        if (npos.lt.2) return
        if (nnew.lt.2) return
        if (nvar.lt.1) return

    !------------- calc. grid-spacing and check uniformity !

        same = (xpos(npos)& 
             -  xpos(  +1)) / (npos-1)

        uniform = .true.
             
        xtol = same * RTOL

        do  ipos = +1 , npos-1, +1
        
            spac = xpos(ipos+1) &
        &        - xpos(ipos+0)
  
            diff = abs(spac - same)
        
            if (diff.gt.xtol) then
            
                uniform = .false.
            
            end if
  
            work% &
        &   cell_spac(ipos) = spac
        
        end do

       !uniform = .false.
    
    !------------- reconstruct FHAT over all cells in XPOS !

        if (.not. uniform) then

    !------------------------------------ variable spacing !
        call rcon1d(npos,nvar,ndof, &
        &           work%cell_spac, &
        &           fdat,bclo,bcup, &
        &           work%cell_func, &
        &           work,opts,tCPU)

        else
        
    !------------------------------------ constant spacing !
        delx(1) = work%cell_spac(1)
        
        call rcon1d(npos,nvar,ndof, &
        &           delx,    &
        &           fdat,bclo,bcup, &
        &           work%cell_func, &
        &           work,opts,tCPU)
        
        end if

    !------------- remap FDAT from XPOS to XNEW using FHAT !

        __TIC__

        select case(opts%cell_meth)
        case(pcm_method)
    !------------------------------------ 1st-order method !
        call imap1d(npos,nnew,nvar, &
        &           ndof,  +1,      &
        &           xpos,xnew,      &
        &           work%cell_func, &
        &           fdat,fnew,xtol)
    
        case(plm_method)
    !------------------------------------ 2nd-order method !
        call imap1d(npos,nnew,nvar, &
        &           ndof,  +2,      &
        &           xpos,xnew,      &
        &           work%cell_func, &
        &           fdat,fnew,xtol)

        case(ppm_method)
    !------------------------------------ 3rd-order method !
        call imap1d(npos,nnew,nvar, &
        &           ndof,  +3,      &
        &           xpos,xnew,      &
        &           work%cell_func, &
        &           fdat,fnew,xtol)

        case(pqm_method)
    !------------------------------------ 5th-order method !
        call imap1d(npos,nnew,nvar, &
        &           ndof,  +5,      &
        &           xpos,xnew,      &
        &           work%cell_func, &
        &           fdat,fnew,xtol)
        
        end select
        
        __TOC__(tCPU,rmap_time)
        
        return
    
    end  subroutine
    
    !------------ IMAP1D: 1-dimensional degree-k remapping !
    
    pure subroutine imap1d(npos,nnew,nvar,ndof, &
        &                  mdof,xpos,xnew,fhat, &
        &                  fdat,fnew,XTOL)

    !
    ! NPOS  no. edges in old grid.
    ! NNEW  no. edges in new grid.
    ! NVAR  no. discrete variables to remap.
    ! NDOF  no. degrees-of-freedom per cell.
    ! MDOF  no. degrees-of-freedom per FHAT.    
    ! XPOS  old grid edge positions. XPOS is a length NPOS
    !       array .
    ! XNEW  new grid edge positions. XNEW is a length NNEW
    !       array .
    ! FHAT  reconstruction over old grid. FHAT has SIZE =
    !       MDOF-by-NVAR-by-NPOS-1 .
    ! FDAT  grid-cell moments on old grid. FDAT has SIZE = 
    !       NDOF-by-NVAR-by-NPOS-1 .
    ! FNEW  grid-cell moments on new grid. FNEW has SIZE = 
    !       NDOF-by-NVAR-by-NNEW-1 .
    ! XTOL  min. grid-cell thickness . 
    !

        implicit none    
    
    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nnew
        integer, intent( in) :: nvar
        integer, intent( in) :: ndof,mdof
        real*8 , intent( in) :: xpos(:)
        real*8 , intent( in) :: xnew(:)
        real*8 , intent( in) :: fhat(:,:,:)
        real*8 , intent( in) :: fdat(:,:,:)        
        real*8 , intent(out) :: fnew(:,:,:)
        real*8 , intent( in) :: XTOL
        
    !------------------------------------------- variables !
        integer :: kpos,ipos,ivar,idof
        integer :: pos0,pos1,vmin,vmax
        real*8  :: xmid,xhat,khat,stmp
        real*8  :: xxlo,xxhi,sslo,sshi,intf
        real*8  :: vvlo(  +1:+5)
        real*8  :: vvhi(  +1:+5)        
        real*8  :: ivec(  +1:+5)
        real*8  :: sdat(  +1:nvar)
        real*8  :: snew(  +1:nvar)
        real*8  :: serr(  +1:nvar)
        integer :: kmin(  +1:nvar)
        integer :: kmax(  +1:nvar)
        
        integer, parameter :: INTB = -1   ! integral basis  

    !------------- remap FDAT from XPOS to XNEW using FHAT !

        kmin = +1 ; kmax = +1
        pos0 = +1 ; pos1 = +1

        do  kpos = +1, nnew-1

    !------ first cell in XPOS overlapping with XNEW(KPOS) !

            pos1 = max(pos1,1)

            do  pos0 = pos1, npos-1

                if (xpos(pos0+1)&
            &       .gt. xnew(kpos+0)) exit

            end do

    !------ final cell in XPOS overlapping with XNEW(KPOS) !
    
            do  pos1 = pos0, npos-1

                if (xpos(pos1+0)&
            &       .ge. xnew(kpos+1)) exit
    
            end do
            
            pos1 = pos1 - 1

    !------------- integrate FHAT across overlapping cells !

            khat = xnew(kpos+1) &
            &    - xnew(kpos+0)
            khat = max (khat , XTOL)

            do  idof = +1,ndof
            do  ivar = +1,nvar
                
                fnew(idof,ivar,kpos) = 0.d0
            
            end do
            end do

            do  ipos = pos0, pos1

    !------------------------------- integration endpoints !
    
                xxlo = max (xpos(ipos+0) , &
            &               xnew(kpos+0))
                xxhi = min (xpos(ipos+1) , &
            &               xnew(kpos+1))

    !------------------------------- local endpoint coords !
    
                xmid = xpos(ipos+1) * .5d0 &
            &        + xpos(ipos+0) * .5d0    
                xhat = xpos(ipos+1) * .5d0 &
            &        - xpos(ipos+0) * .5d0
     
                sslo = &
            &  (xxlo-xmid) / max(xhat,XTOL)
                sshi = &
            &  (xxhi-xmid) / max(xhat,XTOL)

    !------------------------------- integral basis vector !
    
                call bfun1d(INTB,mdof, &
                            sslo,vvlo)
                call bfun1d(INTB,mdof, &
                            sshi,vvhi)
                
                ivec =  vvhi - vvlo

    !--------- integrate FHAT across the overlap XXLO:XXHI !
    
                do  ivar = +1, nvar
    
                intf =  dot_product (  &
            &       ivec(+1:mdof),  &
            &   fhat(+1:mdof,ivar,ipos-0) )

                intf =  intf * xhat
        
    !--------- accumulate integral contributions from IPOS !
    
                fnew(  +1,ivar,kpos) = &
            &   fnew(  +1,ivar,kpos) + intf

                end do

            end do

    !------------------------------- finalise KPOS profile !

            do  ivar = +1, nvar

                fnew(  +1,ivar,kpos) = &
            &   fnew(  +1,ivar,kpos) / khat

    !--------- keep track of MIN/MAX for defect correction !

                vmax =    kmax(ivar)
                vmin =    kmin(ivar)

                if(fnew(1,ivar,kpos) &
            &  .gt.fnew(1,ivar,vmax) ) then
                
                kmax(ivar) =   kpos
            
                else &
            &   if(fnew(1,ivar,kpos) &
            &  .lt.fnew(1,ivar,vmin) ) then
                
                kmin(ivar) =   kpos
            
                end if

            end do

        end do

    !--------- defect corrections: Kahan/Babuska/Neumaier. !

    !   Carefully compute column sums, leading to a defect
    !   wrt. column-wise conservation. Use KBN approach to
    !   account for FP roundoff.

        sdat = 0.d0; serr = 0.d0
        do  ipos = +1, npos-1
        do  ivar = +1, nvar-0
        
    !------------------------------- integrate old profile !

            xhat = xpos(ipos+1) &
        &        - xpos(ipos+0)

            intf = xhat*fdat(1,ivar,ipos)
            
            stmp = sdat(ivar) + intf

            if (abs(sdat(ivar)) &
        &           .ge. abs(intf)) then

            serr(ivar) = &
        &   serr(ivar) + ((sdat(ivar)-stmp)+intf)

            else

            serr(ivar) = &
        &   serr(ivar) + ((intf-stmp)+sdat(ivar))

            end if

            sdat(ivar) = stmp
        
        end do
        end do

        sdat =  sdat + serr

        snew = 0.d0; serr = 0.d0
        do  ipos = +1, nnew-1
        do  ivar = +1, nvar-0
        
    !------------------------------- integrate new profile !

            khat = xnew(ipos+1) &
        &        - xnew(ipos+0)

            intf = khat*fnew(1,ivar,ipos)
            
            stmp = snew(ivar) + intf

            if (abs(snew(ivar)) &
        &           .ge. abs(intf)) then

            serr(ivar) = &
        &   serr(ivar) + ((snew(ivar)-stmp)+intf)

            else

            serr(ivar) = &
        &   serr(ivar) + ((intf-stmp)+snew(ivar))

            end if

            snew(ivar) = stmp
        
        end do
        end do

        snew =  snew + serr
        serr =  sdat - snew

    !--------- defect corrections: nudge away from extrema !

    !   Add a correction to remapped state to impose exact
    !   conservation. Via sign(correction), nudge min/max.
    !   cell means, such that monotonicity is not violated 
    !   near extrema...
    
        do  ivar = +1, nvar-0

            if (serr(ivar) .gt. 0.d0) then

                vmin = kmin(ivar)

                fnew(1,ivar,vmin) = &
        &       fnew(1,ivar,vmin) + &
        &  serr(ivar)/(xnew(vmin+1)-xnew(vmin+0))

            else &
        &   if (serr(ivar) .lt. 0.d0) then

                vmax = kmax(ivar)

                fnew(1,ivar,vmax) = &
        &       fnew(1,ivar,vmax) + &
        &  serr(ivar)/(xnew(vmax+1)-xnew(vmax+0))

            end if

        end do
       
    !------------------------------- new profile now final !

        return
    
    end  subroutine



