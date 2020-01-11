module scheme_mod

use namelist_mod
use mesh_mod

implicit none

private

public full_periodic_boundary_condition
public diffusion_ordinary
public diffusion_flux_limiter

contains

subroutine diffusion_ordinary(w_old, w_tend)
  ! ∂y/∂t = (-1)^((n/2)+1) * a *∂^n y / ∂x^n
  ! a = 1 / dt * (dx/2)^n
  real, intent(in)  :: w_old(:)
  real, intent(out) :: w_tend(:)

  real w_tmp (1-halo: nx+halo),&
       w_tmp2(1-halo: nx+halo),&
       w_new (1-halo: nx+halo)

  integer i, order

  w_tmp(:) = w_old(:)
  w_tmp2(:) = w_old(:)

  do order = 1, diffusion_order / 2
    do i = 1, nx
       w_tend(i) = (w_tmp(i+1) - 2 * w_tmp(i) + w_tmp(i-1)) / (mesh%dx**2)
    end do
    if(order /= diffusion_order / 2) then
       w_tmp(1:nx) = w_tend(1:nx)
    end if
  enddo    

end subroutine diffusion_ordinary

subroutine diffusion_flux_limiter(w_old, w_tend)
  ! xue ming (2000)
  real, intent(in) :: w_old(:)
  real, intent(out) :: w_tend(:)
  real w_tmp(1-halo: nx+halo), &
       w_tmp2(1-halo: nx+halo)
  real full_flux(1-halo: nx+halo)
  real half_grad(nx), half_w_grad(nx)
  integer :: sign, i, order
 
  w_tmp = w_old
  w_tmp2 = w_old

  do order = 1, diffusion_order / 2 -1
    ! 1) Calculate diffusion flux on cell center
    do i = 1, nx
      full_flux(i) = (w_tmp(i+1) - 2 * w_tmp(i) + w_tmp(i-1)) / (mesh%dx**2)
    end do
    call full_periodic_boundary_condition(full_flux)
    if(order /= diffusion_order/2 -1) then
      w_tmp(:) = full_flux(:)
    end if
  end do
  ! 2) Calculate gradient of 2nd/4th order diffusion flux on cell boundary
  do i = 1, nx
     half_grad(i) = (full_flux(i+1) - full_flux(i)) / mesh%dx
  end do 
  
  ! 3) Calculate gradient of orign data on cell boundary
  do i = 1, nx 
     half_w_grad(i) = (w_tmp2(i+1) - w_tmp2(i)) / mesh%dx
  end do 

  ! 4) decide the sign on the cell bounday
  sign = (-1)**(diffusion_order / 2 + 1)
  do i = 1, nx
    if (sign*half_grad(i) * half_w_grad(i) <= 0) then
      half_grad(i) = 0
    end if
  end do 
  call full_periodic_boundary_condition(half_grad)
  do i = 1, nx
     w_tend(i) = (half_grad(i) - half_grad(i-1)) / mesh%dx
  end do

end subroutine diffusion_flux_limiter

subroutine full_periodic_boundary_condition(rho)
  real, intent(inout) :: rho(1-halo:nx+halo)
  integer :: i

  do i = 1, halo
    rho(1-i) = rho(nx-i+1)
    rho(nx+i) = rho(1+i-1)
  end do 
end subroutine full_periodic_boundary_condition

end module scheme_mod