module state_mod

  use namelist_mod
  use mesh_mod
  use scheme_mod

  implicit none

  private

  public state_type
  public state
  public state_init
  public state_final

  type state_type 
    real, allocatable :: rho(:)
  end type state_type

  type(state_type) :: state

contains 

  subroutine state_init()

    integer :: i 
    allocate(state%rho(1-halo: nx+halo))

	  do i = 1, nx
	    if(mesh%x(i) >= 0.28 .and. mesh%x(i) <= 0.72) then 
        state%rho(i) = 0.4 
	    else
        state%rho(i) = -0.4
	    end if 
	  end do 
    call full_periodic_boundary_condition(state%rho)
  end subroutine state_init

  subroutine state_clear(this)
    class(state_type), intent(inout) :: this

    if(allocated(this%rho)) deallocate(this%rho)
  end subroutine state_clear
  
  subroutine state_final()
 
    if(allocated(state%rho))     deallocate(state%rho)
  end subroutine state_final

end module state_mod