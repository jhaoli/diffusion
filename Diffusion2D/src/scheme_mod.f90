module scheme_mod
  
  use namelist_mod
  use mesh_mod
  use boundary_mod
  
  implicit none

  private

  public diffusion_split
  ! public diffusion_direct
  public diffusion_limiter

contains

  subroutine diffusion_split(h_old, h_tend)

    real, intent(in)  :: h_old(1-halo:nx+halo, 1-halo:ny+halo)
    real, intent(out) :: h_tend(1:nx, 1:ny)
 
    real, allocatable :: hx_tend(:,:), &
                         hy_tend(:,:)
                        
    real, allocatable :: hx_tmp(:,:), &
                         hy_tmp(:,:) 
   
    integer :: i, j, order, sign 
 
    allocate(hx_tend(1:nx, 1:ny))
    allocate(hy_tend(1:nx, 1:ny))
    allocate(hx_tmp(1-halo:nx+halo, 1-halo:ny+halo))
    allocate(hy_tmp(1-halo:nx+halo, 1-halo:ny+halo))
 
    hx_tmp = h_old
    hy_tmp = h_old

    do j = 1, ny
      do order = 1, diffusion_order / 2
        do i = 1, nx
          hx_tend(i,j) = (hx_tmp(i+1,j) - 2 * hx_tmp(i,j) + hx_tmp(i-1,j)) / (mesh%dx**2)            
        end do
        if(order /= diffusion_order / 2) then
          hx_tmp(1:nx,j) = hx_tend(1:nx,j)
        end if
      end do
    end do

    do i = 1, nx
      do order = 1, diffusion_order / 2
        do j = 1, ny
          hy_tend(i,j) = (hy_tmp(i,j+1) - 2 * hy_tmp(i,j) + hy_tmp(i,j-1)) / (mesh%dy**2)            
        end do
        if(order /= diffusion_order / 2) then
          hy_tmp(i,1:ny) = hy_tend(i,1:ny)
        end if
      end do
    end do
    h_tend= hx_tend + hy_tend

    deallocate(hx_tmp,hy_tmp, hx_tend, hy_tend)
  end subroutine diffusion_split

  subroutine diffusion_direct(h_old, h_tend)
    
  real, intent(in)  :: h_old(1-halo:nx+halo, 1-halo:ny+halo)
  real, intent(out) :: h_tend(1:nx, 1:ny)
  integer :: i, j, order
  real :: h_tmp(1-halo:nx+halo, 1-halo:ny+halo)

  h_tmp(:,:) = h_old(:,:)
  do order = 1, diffusion_order / 2
    do j = 1, ny
      do i = 1, nx
        h_tend(i,j) = (h_tmp(i+1,j) - 2 * h_tmp(i,j) + h_tmp(i-1,j)) / (mesh%dx**2) +&
                      (h_tmp(i,j+1) - 2 * h_tmp(i,j) + h_tmp(i,j-1)) / (mesh%dy**2)
      end do
    end do 
    if(order /= diffusion_order / 2) then
      h_tmp(1:nx, 1:ny) = h_tend(1:nx, 1:ny)
    end if
  end do
  end subroutine diffusion_direct

  subroutine diffusion_limiter(h_old, h_tend)
  
    real, intent(in)  :: h_old(1-halo:nx+halo, 1-halo:ny+halo)
    real, intent(out) :: h_tend(1:nx, 1:ny)
    real, allocatable :: hx_tend(:,:), &
                         hy_tend(:,:)
                        
    real, allocatable :: hx_tmp(:,:), &
                         hy_tmp(:,:) 
    real, allocatable :: half_grad_x(:,:), &
                         half_grad_y(:,:), &
                         half_h_grad_x(:,:), &
                         half_h_grad_y(:,:)
    integer :: i, j, order, sign 
    allocate(hx_tend(1-halo:nx+halo, 1-halo:ny+halo))
    allocate(hy_tend(1-halo:nx+halo, 1-halo:ny+halo))
    allocate(hx_tmp (1-halo:nx+halo, 1-halo:ny+halo))
    allocate(hy_tmp (1-halo:nx+halo, 1-halo:ny+halo))
    allocate(half_grad_x  (0:nx, 1:ny))
    allocate(half_grad_y  (1:nx, 0:ny))
    allocate(half_h_grad_x(0:nx, 1:ny))
    allocate(half_h_grad_y(1:nx, 0:ny))
  
      ! 1) Calculate even order laplacian
    hx_tmp(:,:) = h_old(:,:)
    do j = 1, ny
      do order = 1, diffusion_order / 2 -1
        do i = 0, nx+1
           hx_tend(i,j) = (hx_tmp(i+1,j) - 2 * hx_tmp(i,j) + hx_tmp(i-1,j)) / (mesh%dx**2)            
        end do
        if(order /= diffusion_order / 2 - 1) then
           hx_tmp(:,j) = hx_tend(:,j)
        end if
      end do
    end do
  
    hy_tmp(:,:) = h_old(:,:)
    do i = 1, nx
      do order = 1, diffusion_order / 2 -1
        do j = 0, ny+1
           hy_tend(i,j) = (hy_tmp(i,j+1) - 2 * hy_tmp(i,j) + hy_tmp(i,j-1)) / (mesh%dy**2)            
        end do
        if(order /= diffusion_order / 2 - 1) then
           hy_tmp(i,:) = hy_tend(i,:)
        end if
      end do
    end do 
  
    ! 2) Calculate gradient of laplcain above 
    do j = 1, ny
      do i = 0, nx
        half_grad_x(i,j) = (hx_tend(i+1,j) - hx_tend(i,j)) / mesh%dx
      end do 
    end do 
    do i = 1, nx
      do j = 0, ny
        half_grad_y(i,j) = (hy_tend(i,j+1) - hy_tend(i,j)) / mesh%dy
      end do 
    end do
     
    ! 3) Calculate gradient of the variable be diffused 
    do j = 1, ny
      do i = 0, nx
        half_h_grad_x(i,j) = (h_old(i+1,j) - h_old(i,j)) / mesh%dx
      end do 
    end do 
    do i = 1, nx
      do j = 0, ny
        half_h_grad_y(i,j) = (h_old(i,j+1) - h_old(i,j)) / mesh%dy
      end do 
    end do
  
    !4) decide the sign of diffusion flux to be zero or not
    sign = (-1)**(diffusion_order / 2 + 1)
    do j = 1, ny
      do i = 0, nx
        if (sign * half_grad_x(i,j) * half_h_grad_x(i,j) <= 0) then
          half_grad_x(i,j) = 0.0
        end if 
      end do 
    end do 
    do i = 1, nx
      do j = 0, ny
        if (sign * half_grad_y(i,j) * half_h_grad_y(i,j) <= 0) then
          half_grad_y(i,j) = 0.0
        end if
      end do 
    end do
  
    !5) Calculate the time tendency lastly.
    do j = 1, ny
      do i = 1, nx
        hx_tend(i,j) = (half_grad_x(i,j) - half_grad_x(i-1,j)) / mesh%dx
      end do
    end do
    do i = 1, nx
      do j = 1, ny
        hy_tend(i,j) = (half_grad_y(i,j) - half_grad_y(i,j-1)) / mesh%dy
      end do 
    end do
    
    do j = 1, ny
      do i = 1, nx
        h_tend(i,j) = hx_tend(i,j) + hy_tend(i,j)
      end do 
    end do 
    
    deallocate(hx_tend, hy_tend, hx_tmp, hy_tmp)
    deallocate(half_grad_x, half_grad_y, half_h_grad_x, half_h_grad_y)
  end subroutine diffusion_limiter

end module scheme_mod


!! directly apply the 4th order diffusion
! do j = 2, num_points-1
!    do i = 2, num_points-1
!       h_tend(i,j) = (hx_tmp(i+2,j) - 4 * hx_tmp(i+1,j) + 6 * hx_tmp(i,j) - 4 * hx_tmp(i-1,j) + hx_tmp(i-2,j)) / (dx**4) +&
!                     (hy_tmp(i,j+2) - 4 * hy_tmp(i,j+1) + 6 * hy_tmp(i,j) - 4 * hy_tmp(i,j-1) + hy_tmp(i,j-2)) / (dy**4)
!    end do
! end do
!! end 
! h_tend(0:1,:) = 0
! h_tend(num_points:num_points+1,:) = 0
! h_tend(:,0:1) = 0
! h_tend(:, num_points:num_points+1) = 0