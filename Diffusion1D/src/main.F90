program main

  use namelist_mod
  use mesh_mod
  use state_mod
  use scheme_mod
  use history_mod

  implicit none

  integer :: step
  character(256) :: namelist_file_path

  interface
    subroutine diffusion_interface(rho)
      real, intent(inout) :: rho(:) 
    end subroutine
  end interface 
  procedure(diffusion_interface), pointer :: diffuse

!!! begin
  call get_command_argument(1, namelist_file_path)
  if (namelist_file_path =='') then
     print*, 'You should give a namelist file path!'
     stop
  end if 
  call parse_namelist(namelist_file_path)

  select case(diffusion_method)
  case ('ordinary')
    diffuse => diffusion_ordinary
  case ('flux_limiter')
    diffuse => diffusion_flux_limiter
  case default
    write(*,*) 'Unknown diffusion method' // trim(diffusion_method) // '!'
  end select

  call mesh_init()
  call state_init()
  call output(0, state)

  do step = 1, num_of_time_step
    call diffuse(state%rho)
    if (mod(step, freq_output) == 0) then
      call output(step, state)
    endif
  enddo

   ! Free spaces
  call mesh_final()
  call state_final()

end program main
