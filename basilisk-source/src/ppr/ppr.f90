module my_module

  use ppr_1d
  use :: iso_c_binding
  implicit none

contains

  subroutine my_remap(npos,nnew,nvar,ndof,xpos,xnew,fdat,fnew,edge_meth,cell_meth, cell_lim) bind(C,name='my_remap')

    integer, intent( in) :: npos,nnew
    integer, intent( in) :: nvar,ndof
    integer, intent( in) :: edge_meth,cell_meth,cell_lim
    real(c_double) , intent( in), dimension(npos) :: xpos
    real(c_double) , intent( in), dimension(nnew) :: xnew
    real(c_double) , intent( in), dimension(ndof,nvar,(npos-1)) :: fdat
    real(c_double) , intent(out), dimension(ndof,nvar,(nnew-1)) :: fnew

    type(rmap_work) :: work  
    type(rmap_opts) :: opts
    type(rcon_ends) :: bc_l(nvar)
    type(rcon_ends) :: bc_r(nvar)

    ! !------------------------------ specify method options !

    opts%edge_meth = edge_meth
    opts%cell_meth = cell_meth
    opts%cell_lims = cell_lim

    ! !------------------------------ set BC.'s at endpoints !

    bc_l%bcopt = bcon_loose         ! "loose" = extrapolate
    bc_r%bcopt = bcon_loose

    ! !------------------------------ init. method workspace !

    call work%init(npos,nvar,opts)

    call rmap1d(npos,nnew,nvar,ndof, &
          &           xpos,xnew,fdat,fnew, &
          &           bc_l,bc_r,work,opts)

    call work%free()
  end subroutine my_remap
  
end module my_module
