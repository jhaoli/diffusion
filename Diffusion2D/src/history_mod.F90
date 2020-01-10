module history_mod

use namelist_mod
use netcdf

implicit none 

private 

public output

contains

  subroutine output(x, y, time_step, rho)
    real, intent(in) :: x(nx), y(ny)
    real, intent(in) :: rho(1:nx, 1:ny)
    integer, intent(in) :: time_step
    integer, parameter :: ndims = 2
  
    character(30) file_name
    integer file_id, time_dim_id, time_var_id, x_dim_id, y_dim_id,&
            x_var_id, y_var_id, rho_var_id, ierr
    integer :: dimids(ndims)
  
    write(file_name, "('diffusion.', I3.3, '.', I4.4, '.nc')") nx, time_step
    
    call check(nf90_create(file_name, nf90_clobber, file_id))
    
    call check(nf90_def_dim(file_id, 'time', nf90_unlimited, time_dim_id))
    call check(nf90_def_var(file_id, 'time', nf90_int, [time_dim_id], time_var_id))
    
    call check(nf90_def_dim(file_id, 'x', nx, x_dim_id))
    call check(nf90_def_var(file_id, 'x', nf90_real, [x_dim_id], x_var_id))
  
    call check(nf90_def_dim(file_id, 'y', nx, y_dim_id))
    call check(nf90_def_var(file_id, 'y', nf90_real, [y_dim_id], y_var_id))
  
    dimids = (/x_dim_id, y_dim_id/)
    call check(nf90_def_var(file_id, 'rho', nf90_real, dimids, rho_var_id))
  
    call check(nf90_enddef(file_id))
  
    call check(nf90_put_var(file_id, time_var_id, time_step))

    call check(nf90_put_var(file_id, x_var_id, x))
    call check(nf90_put_var(file_id, y_var_id, y))

    call check(nf90_put_var(file_id, rho_var_id, rho(1:ny,1:nx)))
 
    call check(nf90_close(file_id))
  end subroutine output

  subroutine check(status)
    integer, intent(in) :: status
  
    if (status/=nf90_noerr) then
      write(*, *) trim(nf90_strerror(status))
      stop "Stopped"
    end if 
  end subroutine check

end module history_mod

