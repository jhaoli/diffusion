program main

  use namelist_mod
  use mesh_mod
  use state_mod
  use scheme_mod
  use history_mod

  implicit none
  
  real, allocatable :: w_new(:)
  real, allocatable :: w_tend(:)
  character(256) :: namelist_file_path
  real :: diffusion_coef
  integer :: i, step, order, sign

  ! interface
  !   subroutine diffusion_interface(rho, rho_tend)
  !     real, intent(in) :: rho(:) 
  !     real, intent(out) :: rho_tend(:)
  !   end subroutine
  ! end interface 
  ! procedure(diffusion_interface), pointer :: diffuse

!!! begin
  call get_command_argument(1, namelist_file_path)
  if (namelist_file_path =='') then
     print*, 'You should give a namelist file path!'
     stop
  end if 
  call parse_namelist(namelist_file_path)

  ! select case(diffusion_method)
  ! case ('ordinary')
  !   diffuse => diffusion_ordinary
  ! case ('flux_limiter')
  !   diffuse => diffusion_flux_limiter
  ! case default
  !   write(*,*) 'Unknown diffusion method' // trim(diffusion_method) // '!'
  ! end select

  call mesh_init()
  call state_init()
  call output(0, state)
  
  allocate(w_tend(1:nx))
  allocate(w_new(1-halo:nx+halo))

  do step = 1, num_of_time_step
    select case(diffusion_method)
    case ('ordinary')
      call diffusion_ordinary(state%rho, w_tend)
    case ('flux_limiter')
      call diffusion_flux_limiter(state%rho, w_tend)
    case default
      write(*,*) 'Unknown diffusion method' // trim(diffusion_method) //'!'
    end select
    ! call diffuse(state%rho, w_tend)
    diffusion_coef = 1.0 / time_step_size * (mesh%dx/2.0)**diffusion_order
    sign = (-1)**(diffusion_order / 2 + 1)
    do i = 1, nx+1
       w_new(i) = state%rho(i) + sign * time_step_size * diffusion_coef * w_tend(i)
    end do
    call full_periodic_boundary_condition(w_new)
    state%rho(:) = w_new(:)
    if (mod(step, freq_output) == 0) then
      call output(step, state)
    endif
  enddo

   ! Free spaces
  call mesh_final()
  call state_final()
  deallocate(w_tend)
  deallocate(w_new)
end program main
