module state_mod

  use namelist_mod
  use mesh_mod
  use boundary_mod

  implicit none

  private
  public state_type
  public state_init
  public state

  type state_type
    real, allocatable :: h(:,:)
  end type 

  type(state_type) :: state
contains
  
  subroutine state_init()
    integer :: i, j 
    allocate(state%h(1-halo:nx+halo, 1-halo:ny+halo))

    state%h = -0.4
    do j = 1, ny
      do i = 1, nx
        if (mesh%x(i)>= 0.3 .and. mesh%x(i) <=0.7 .and. mesh%y(j) >=0.3 .and. mesh%y(j)<=0.7) then
          state%h(i,j) = 0.4
        end if 
      end do
    end do
    call full_periodic_boundary_condition(state%h)

  end subroutine state_init

end module state_mod