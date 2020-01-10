module namelist_mod

  implicit none 

  real :: time_step_size = 1
  real :: xmin = 0.0 
  real :: xmax = 1.0
  integer :: num_of_time_step = 100 
  integer :: freq_output = 10
  integer :: halo = 4
  integer :: nx = 50
  integer :: diffusion_order = 8
  character(20) :: diffusion_method = 'flux_limiter' !  'split' or 'direct' or 'limiter'

  namelist /params/     &
      time_step_size   ,&
      num_of_time_step ,& 
      freq_output      ,& ! output frequency(every so many steps)
      xmin             ,&
      xmax             ,&
      halo             ,&
      nx               ,&
      diffusion_method ,&
      diffusion_order   

contains

  subroutine parse_namelist(file_path)

    character(*), intent(in) :: file_path

    open(10, file=file_path, status='old')
    read(10, nml=params)
    close(10)

  end subroutine parse_namelist

end module namelist_mod