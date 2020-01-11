program main
  
  use namelist_mod
  use mesh_mod
  use state_mod
  use scheme_mod
  use history_mod
  use boundary_mod

  implicit none

  real, allocatable :: h_new(:,:)
  real, allocatable :: h_tend(:,:)
  real :: diffusion_coef
  integer :: i, j, step, order, sign
  
  character(256) :: namelist_file_path
  ! interface
  !   subroutine diffusion_interface(h, rho_tend)
  !     real, intent(in) :: h(:,:) 
  !     real, intent(out) :: rho_tend(:,:)
  !   end subroutine
  ! end interface 
  ! procedure(diffusion_interface), pointer :: diffuse 

  call get_command_argument(1, namelist_file_path)
  if (namelist_file_path =='') then
    print*, 'You should give a namelist file path!'
    stop
  end if 
  call parse_namelist(namelist_file_path)
  
  ! select case(diffusion_method)
  ! case ('split')
  !   diffuse => diffusion_split
  ! case ('limiter')
  !   diffuse => diffusion_limiter
  ! case default
  !   write(*,*) 'Unknown diffusion method' // trim(diffusion_method) // '!'
  ! end select

  allocate(h_new(1-halo:nx+halo, 1-halo:ny+halo))
  allocate(h_tend(1:nx, 1:ny))
 
  call mesh_init()
  call state_init()
  call output(0, state%h(1:nx, 1:ny))

  do step = 1, num_of_time_step
    select case(diffusion_method)
    case ('split')
      call diffusion_split(state%h, h_tend)
    case ('limiter')
      call diffusion_limiter(state%h, h_tend)
    case default
      write(*,*) 'Unknown diffusion method' // trim(diffusion_method) // '!'
    end select
    ! call diffusion_split(state%h, h_tend)
    diffusion_coef = 1.0 / time_step_size * (mesh%dx/2.0)**diffusion_order
    sign = (-1)**(diffusion_order / 2 + 1)
    
    do j = 1, ny
      do i = 1, nx
      h_new(i,j) = state%h(i,j) + sign * time_step_size * diffusion_coef * h_tend(i,j)
      end do
    end do
    call full_periodic_boundary_condition(h_new)
    state%h(:,:) = h_new(:,:)     
    if (mod(step, freq_output) == 0) then
      call output(step, h_new(1:nx, 1:ny))
    endif
  enddo
  deallocate(h_new, h_tend)

end program main
