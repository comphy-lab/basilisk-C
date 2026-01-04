
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

    module ppr_1d

    !    
    ! PPR-1D.f90: 1-d piecewise polynomial reconstructions.
    !
    ! Darren Engwirda 
    ! 31-Mar-2019
    ! de2363 [at] columbia [dot] edu
    !
    !

    implicit none

    !------------------------------------ compile-time def !

!   define __PPR_WARNMAT__
!   define __PPR_DOTIMER__

#   ifdef  __PPR_DOTIMER__

#       define __TIC__                      \
        call system_clock   (ttic,rate)

#       define __TOC__(time, mark)          \
        call system_clock   (ttoc,rate) ;   \
        if ( present(time)) \
        time%mark=time%mark+(ttoc-ttic)
    
#   else

#       define __TIC__
#       define __TOC__(time, mark)

#   endif

    !------------------------------------ method selection !
        
    integer, parameter :: p1e_method = +100
    integer, parameter :: p3e_method = +101
    integer, parameter :: p5e_method = +102

    integer, parameter :: pcm_method = +200
    integer, parameter :: plm_method = +201
    integer, parameter :: ppm_method = +202
    integer, parameter :: pqm_method = +203

    integer, parameter :: null_limit = +300
    integer, parameter :: mono_limit = +301
    integer, parameter :: weno_limit = +302

    integer, parameter :: bcon_loose = +400
    integer, parameter :: bcon_value = +401
    integer, parameter :: bcon_slope = +402
 
    type rmap_tics
    !------------------------------- tCPU timer for RCON1D !
        integer :: rmap_time
        integer :: edge_time
        integer :: cell_time
        integer :: oscl_time
    end type rmap_tics

    type rcon_opts
    !------------------------------- parameters for RCON1D !
        integer :: edge_meth
        integer :: cell_meth
        integer :: cell_lims
        integer :: wall_lims
    end type rcon_opts

    type rcon_ends
    !------------------------------- end-conditions struct !
        integer :: bcopt
        real*8  :: value
        real*8  :: slope
    end type rcon_ends

    type rcon_work
    !------------------------------- work-space for RCON1D !
        real*8, allocatable  :: edge_func(:,:)
        real*8, allocatable  :: edge_dfdx(:,:)
        real*8, allocatable  :: cell_oscl(:,:,:)
    contains
        procedure :: init => init_rcon_work
        procedure :: free => free_rcon_work
    end type rcon_work

    type, extends(rcon_opts) :: rmap_opts
    !------------------------------- parameters for RMAP1D !
    end type rmap_opts

    type, extends(rcon_work) :: rmap_work
    !------------------------------- work-space for RMAP1D !
        real*8, allocatable  :: cell_spac(:)
        real*8, allocatable  :: cell_func(:,:,:)
    contains
        procedure :: init => init_rmap_work
        procedure :: free => free_rmap_work
    end type rmap_work

    contains
 
    !------------------------------------------------------!
    ! INIT-RCON-WORK: init. work-space for RCON1D.         !
    !------------------------------------------------------!
    
    subroutine init_rcon_work(this,npos,nvar,opts)

    !
    ! THIS  work-space structure for RCON1D .
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! OPTS  parameters structure for RCON1D .
    !

        implicit none

    !------------------------------------------- arguments !
        class(rcon_work) , intent(inout) :: this
        integer, intent(in):: npos
        integer, intent(in):: nvar
        class(rcon_opts) , optional      :: opts

    !------------------------------------------- variables !
        integer :: okay

        allocate(this% &
        &   edge_func(  nvar,npos), &
        &        this% &
        &   edge_dfdx(  nvar,npos), &
        &        this% &
        &   cell_oscl(2,nvar,npos), &
        &   stat=okay)
        
    end subroutine

    !------------------------------------------------------!
    ! INIT-RMAP-WORK: init. work-space for RMAP1D.         !
    !------------------------------------------------------!
    
    subroutine init_rmap_work(this,npos,nvar,opts)

    !
    ! THIS  work-space structure for RMAP1D .
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! OPTS  parameters structure for RMAP1D .
    !

        implicit none

    !------------------------------------------- arguments !
        class(rmap_work) , intent(inout) :: this
        integer, intent(in) :: npos
        integer, intent(in) :: nvar
        class(rcon_opts) , optional      :: opts

    !------------------------------------------- variables !
        integer :: okay,ndof

        ndof = ndof1d(opts%cell_meth)

        allocate(this% &
        &   edge_func(  nvar,npos), &
        &        this% &
        &   edge_dfdx(  nvar,npos), &
        &        this% &
        &   cell_oscl(2,nvar,npos), &
        &        this% &
        &   cell_spac(       npos), &
        &        this% &
        &   cell_func(ndof,nvar,npos) , &
        &   stat=okay)
 
    end subroutine

    !------------------------------------------------------!
    ! FREE-RCON-WORK: free work-space for RCON1D .         !
    !------------------------------------------------------!
    
    subroutine free_rcon_work(this)

        implicit none

    !------------------------------------------- arguments !
        class(rcon_work), intent(inout) :: this

        deallocate(this%edge_func, &
        &          this%edge_dfdx, &
        &          this%cell_oscl)
     
    end subroutine

    !------------------------------------------------------!
    ! FREE-RMAP-WORK: free work-space for RMAP1D .         !
    !------------------------------------------------------!
    
    subroutine free_rmap_work(this)

        implicit none

    !------------------------------------------- arguments !
        class(rmap_work), intent(inout) :: this


        deallocate(this%edge_func, &
        &          this%edge_dfdx, &
        &          this%cell_oscl, &
        &          this%cell_func, &
        &          this%cell_spac)

    end subroutine

    !------------------------------------------------------!
    ! NDOF1D : no. degrees-of-freedom per polynomial .     !
    !------------------------------------------------------!
    
    pure function ndof1d(meth) result(rdof)

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: meth

    !------------------------------------------- variables !
        integer  :: rdof

        select case(meth)
    !-------------------------------- edge reconstructions !
        case (p1e_method)
            rdof = +2
        case (p3e_method)
            rdof = +4
        case (p5e_method)
            rdof = +6
    !-------------------------------- cell reconstructions !
        case (pcm_method)
            rdof = +1
        case (plm_method)
            rdof = +2
        case (ppm_method)
            rdof = +3
        case (pqm_method)
            rdof = +5
        
        case default 
            rdof = +0

        end select
        
    end function  ndof1d

    !------------------------------------------------------!
    ! BFUN1D : one-dimensional poly. basis-functions .     !
    !------------------------------------------------------!
        
#       include "bfun1d.f90"

    !------------------------------------------------------!
    ! UTIL1D : one-dimensional grid manip. utilities .     !
    !------------------------------------------------------!
    
#       include "util1d.f90"

    !------------------------------------------------------!
    ! WENO1D : "essentially" non-oscillatory limiter .     !
    !------------------------------------------------------!

#       include "weno1d.f90"

#       include "oscl1d.f90"

    !------------------------------------------------------!
    ! RCON1D : one-dimensional poly. reconstructions .     !
    !------------------------------------------------------!

#       include "rcon1d.f90"

#       include "inv.f90"

#       include "pbc.f90"
#       include "p1e.f90"
#       include "p3e.f90"
#       include "p5e.f90"

#       include "root1d.f90"
    
#       include "pcm.f90"
#       include "plm.f90"
#       include "ppm.f90"
#       include "pqm.f90"
    
    !------------------------------------------------------!
    ! RMAP1D : one-dimensional conservative "re-map" .     !
    !------------------------------------------------------!
       
#       include "rmap1d.f90"

    !------------------------------------------------------!
    ! FFSL1D : one-dimensional FFSL scalar transport .     !
    !------------------------------------------------------!
    
#       include "ffsl1d.f90"


    !------------------------------------------ end ppr_1d !
      
#       undef   __TIC__
#       undef   __TOC__
        
    end module
    
    
    
