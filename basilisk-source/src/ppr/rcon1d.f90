
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
    ! RCON1D.f90: conservative, polynomial reconstructions.
    !
    ! Darren Engwirda 
    ! 07-Sep-2016
    ! de2363 [at] columbia [dot] edu
    !
    !
    
    subroutine rcon1d(npos,nvar,ndof,delx,fdat, &
        &             bclo,bchi,fhat,work,opts, &
        &             tCPU)

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
    ! FHAT  grid-cell re-con. array. FHAT is an array with
    !       SIZE = MDOF-by-NVAR-by-NPOS-1 .
    ! WORK  method work-space. See RCON-WORK for details .
    ! OPTS  method parameters. See RCON-OPTS for details .
    ! TCPU  method tcpu-timer.
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        class(rcon_work), intent(inout):: work
        class(rcon_opts), intent(in)   :: opts
        real*8 , intent( in) :: delx(:)
        real*8 , intent(out) :: fhat(:,:,:)
        real*8 , intent( in) :: fdat(:,:,:)
        type (rcon_ends), intent(in) :: bclo(:)
        type (rcon_ends), intent(in) :: bchi(:)
        type (rmap_tics), &
        &   intent(inout) , optional :: tCPU

    !------------------------------------------- variables !
        integer :: halo,ipos
        real*8  :: dmin,dmid

#       ifdef __PPR_TIMER__
        integer(kind=8) :: ttic,ttoc,rate
#       endif

        if (ndof.lt.1) return
        if (npos.lt.2) return
        if (nvar.lt.1) return
     
    !-------------------------- compute min grid-tolerance !
     
        dmid = delx(1)
     
        if (size(delx).gt.+1) then
        
            do  ipos = 2, npos-1
                dmid = &
            &   dmid + delx (ipos)
            end do
        
            dmid = dmid /(npos-1)
        
        end if

        dmin = +1.0d-14 * dmid
        
    !-------------------------- compute edge values/slopes !

        __TIC__

        if ( (opts%cell_meth.eq.ppm_method) &
    &  .or.  (opts%cell_meth.eq.pqm_method) ) then

        select case (opts%edge_meth)
            case(p1e_method)
    !------------------------------------ 2nd-order method !
            halo = +1
            call p1e(npos,nvar,ndof, &
            &        delx,fdat,      &
            &        bclo,bchi,      &
            &        work%edge_func, &
            &        work%edge_dfdx, &
            &        opts,dmin)

            case(p3e_method)
    !------------------------------------ 4th-order method !           
            halo = +2
            call p3e(npos,nvar,ndof, &
            &        delx,fdat,      &
            &        bclo,bchi,      &
            &        work%edge_func, &
            &        work%edge_dfdx, &
            &        opts,dmin)

            case(p5e_method)
    !------------------------------------ 6th-order method !           
            halo = +3
            call p5e(npos,nvar,ndof, &
            &        delx,fdat,      &
            &        bclo,bchi,      &
            &        work%edge_func, &
            &        work%edge_dfdx, &
            &        opts,dmin)

        end select

        end if
        
        __TOC__(tCPU,edge_time)

    !-------------------------- compute oscil. derivatives !

        __TIC__

        if (opts%cell_lims.eq.weno_limit) then
    
            call oscli(npos,nvar,ndof, &
            &          delx,fdat, &
            &          work%cell_oscl, &
            &          dmin)
     
        end if
    
        __TOC__(time,oscl_time)

    !-------------------------- compute grid-cell profiles !

        __TIC__

        select case (opts%cell_meth)
            case(pcm_method)
    !------------------------------------ 1st-order method !
            call pcm(npos,nvar,ndof, &
            &        fdat,fhat)

            case(plm_method)
    !------------------------------------ 2nd-order method !
            call plm(npos,nvar,ndof, &
            &        delx,fdat,fhat, &
            &        dmin,&
            &        opts%cell_lims)

            case(ppm_method)
    !------------------------------------ 3rd-order method !
            call ppm(npos,nvar,ndof, &
            &        delx,fdat,fhat, &
            &        work%edge_func, &
            &        work%cell_oscl, &
            &        dmin,&
            &        opts%cell_lims, &
            &        opts%wall_lims, &
            &        halo )

            case(pqm_method)
    !------------------------------------ 5th-order method !
            call pqm(npos,nvar,ndof, &
            &        delx,fdat,fhat, &
            &        work%edge_func, &
            &        work%edge_dfdx, &
            &        work%cell_oscl, &
            &        dmin,&
            &        opts%cell_lims, &
            &        opts%wall_lims, &
            &        halo )

        end select

        __TOC__(tCPU,cell_time)

    end subroutine
    
    
    
