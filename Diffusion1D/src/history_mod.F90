module history_mod

use namelist_mod
use state_mod
use mesh_mod
use netcdf

implicit none 

private 

public output

contains

  subroutine output(time_step, state)

    type(state_type), intent(in) :: state
    integer, intent(in) :: time_step
    character(30):: file_name
    integer:: file_id, time_dim_id, time_var_id, x_dim_id,&
            x_var_id, rho_var_id, ierr
    integer::  status

    write(file_name, "('diffusion.', I3.3, '.', I4.4, '.nc')") nx, time_step
    
    status = nf90_create(file_name, nf90_clobber, file_id)
    call check(status, 'open')
    
    status = nf90_def_dim(file_id, 'time', nf90_unlimited, time_dim_id)
    call check(status, 'define dimension of time')
    status = nf90_def_var(file_id, 'time', nf90_int, [time_dim_id], time_var_id)
    call check(status, 'define variable of time')

    
    status = nf90_def_dim(file_id, 'x', nx, x_dim_id)
    call check(status, 'define dimension of x axis')
    status = nf90_def_var(file_id, 'x', nf90_real, [x_dim_id], x_var_id)
    call check(status, 'define variable of x')
  
    status = nf90_def_var(file_id, 'rho', nf90_float, [x_dim_id, time_dim_id], rho_var_id)
    call check(status, 'define variable of rho')
  
    status = nf90_enddef(file_id)
    call check(status,'end define file')
  
    status = nf90_put_var(file_id, time_var_id, time_step)
    call check(status, 'put variable of time')

    status = nf90_put_var(file_id, x_var_id, mesh%x)
    call check(status, 'put variable of x')

    status = nf90_put_var(file_id, rho_var_id, state%rho(1:nx))
    call check(status, 'put varibale of rho')
 
    status = nf90_close(file_id)
    call check(status, 'close')
  end subroutine output

  subroutine check(status, operation)
    
    integer, intent(in) :: status
    character(*), intent(in) :: operation
  
    if (status == nf90_noerr) return
    write(*, *) "Error encountered during ", operation
    write(*, *) trim(nf90_strerror(status))
    stop 1
 
  end subroutine check

end module history_mod

