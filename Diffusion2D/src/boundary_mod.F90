module boundary_mod
  use namelist_mod

  implicit none

  private
  public full_periodic_boundary_condition

contains

  subroutine full_periodic_boundary_condition(f)

    real, intent(inout) :: f(1-halo:nx+halo,1-halo:ny+halo)
    integer :: i, j

    do j = 1, ny
      do i = 1, halo
        f(1-i , j) = f(nx-i+1, j)
        f(nx+i, j) = f(1+i-1 , j)
      end do 
    end do 
    do i = 1, nx
      do j = 1, halo
        f(i, 1-j ) = f(i, ny-j+1)
        f(i, ny+j) = f(i, 1+j-1 )
      end do
    end do
  end subroutine full_periodic_boundary_condition
end module