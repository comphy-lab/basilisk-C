subroutine cvmix_redirect_stdout()
  open(unit=6, file='/dev/stderr', form='formatted')
end subroutine cvmix_redirect_stdout

subroutine allocate1d(len, mem)
  integer, intent(in) :: len
  real(8), dimension(:), pointer, intent(out) :: mem
  allocate(mem(len))
end subroutine allocate1d

subroutine deallocate1d(mem)
  real(8), dimension(:), pointer, intent(inout) :: mem
  deallocate(mem)
end subroutine deallocate1d

subroutine allocate2d(len1, len2, mem)
  integer, intent(in) :: len1, len2
  real(8), dimension(:,:), pointer, intent(out) :: mem
  allocate(mem(len1,len2))
  do i=1,len1
     do j=1,len2
        mem(i,j) = i
     end do
  end do
end subroutine allocate2d

subroutine deallocate2d(mem)
  real(8), dimension(:,:), pointer, intent(inout) :: mem
  deallocate(mem)
end subroutine deallocate2d

subroutine sizeofall()
  use cvmix_kinds_and_types, only : cvmix_r8, cvmix_data_type

!  real(cvmix_r8), dimension(:), pointer :: a
  type(cvmix_data_type) :: b

!  print *, storage_size(a) ! bits
  print *, sizeof(b)       ! bytes
!  print *, storage_size(b) ! bits

  print *, loc(b%max_nlev) - loc(b)
  print *, loc(b%oceandepth) - loc(b)
  print *, loc(b%SimmonsCoeff) - loc(b)
  print *, loc(b%zw_iface) - loc(b)
  print *, loc(b%zt_cntr) - loc(b)
end subroutine sizeofall

