module diffusion_scheme

private

public diffusion_split
public diffusion_direct
public diffusion_limiter
public full_periodic_boundary_condition

contains

subroutine diffusion_split(h_old, h_tend, num_cells, diffusion_order, dx, dy, halo)
  implicit none
   integer, intent(in) :: num_cells, diffusion_order, halo
   real, intent(in)  :: h_old(1-halo:num_cells+halo, 1-halo:num_cells+halo), dx, dy
   real, intent(out) :: h_tend(1:num_cells, 1:num_cells)

   real, allocatable :: hx_tend(:,:), &
                        hy_tend(:,:)
                       
   real, allocatable :: hx_tmp(:,:), &
                        hy_tmp(:,:) 
  
   integer :: i, j, order, sign 

   allocate(hx_tend(1:num_cells, 1:num_cells))
   allocate(hy_tend(1:num_cells, 1:num_cells))
   allocate(hx_tmp(1-halo:num_cells+halo, 1-halo:num_cells+halo))
   allocate(hy_tmp(1-halo:num_cells+halo, 1-halo:num_cells+halo))
 
    hx_tmp(:,:) = h_old(:,:)
    hy_tmp(:,:) = h_old(:,:)

    do j = 1, num_cells
      do order = 1, diffusion_order / 2
        do i = 1, num_cells
          hx_tend(i,j) = (hx_tmp(i+1,j) - 2 * hx_tmp(i,j) + hx_tmp(i-1,j)) / (dx**2)            
        end do
        if(order /= diffusion_order / 2) then
          hx_tmp(1:num_cells,j) = hx_tend(1:num_cells,j)
        end if
      end do
    end do

    do i = 1, num_cells
      do order = 1, diffusion_order / 2
        do j = 1, num_cells
          hy_tend(i,j) = (hy_tmp(i,j+1) - 2 * hy_tmp(i,j) + hy_tmp(i,j-1)) / (dy**2)            
        end do
        if(order /= diffusion_order / 2) then
          hy_tmp(i,1:num_cells) = hy_tend(i,1:num_cells)
        end if
      end do
    end do
    h_tend= hx_tend + hy_tend
    return
    deallocate(hx_tmp,hy_tmp, hx_tend, hy_tend)
end subroutine diffusion_split

subroutine diffusion_direct(h_old, h_tend, num_cells, diffusion_order, dx, dy, halo)
    implicit none
   real, intent(in)  :: h_old(1-halo:num_cells+halo, 1-halo:num_cells+halo), dx, dy
   real, intent(out) :: h_tend(1:num_cells, 1:num_cells)
   integer, intent(in) :: num_cells, diffusion_order, halo
   integer :: i, j, order
   real :: h_tmp(1-halo:num_cells+halo, 1-halo:num_cells+halo)

   h_tmp(:,:) = h_old(:,:)
   do order = 1, diffusion_order / 2
      do j = 1, num_cells
         do i = 1, num_cells
            h_tend(i,j) = (h_tmp(i+1,j) - 2 * h_tmp(i,j) + h_tmp(i-1,j)) / (dx**2) +&
                          (h_tmp(i,j+1) - 2 * h_tmp(i,j) + h_tmp(i,j-1)) / (dy**2)
         end do
      end do 
      if(order /= diffusion_order / 2) then
        h_tmp(1:num_cells, 1:num_cells) = h_tend(1:num_cells, 1:num_cells)
      end if
    end do
    return
end subroutine diffusion_direct

subroutine diffusion_limiter(h_old, h_tend, num_cells, diffusion_order, dx, dy, halo)
   implicit none
   real, intent(in)  :: h_old(1-halo:num_cells+halo, 1-halo:num_cells+halo), dx, dy
   real, intent(out) :: h_tend(1:num_cells, 1:num_cells)
   integer, intent(in) :: num_cells, diffusion_order, halo 

   real, allocatable :: hx_tend(:,:), &
                        hy_tend(:,:)
                       
   real, allocatable :: hx_tmp(:,:), &
                        hy_tmp(:,:) 
   real, allocatable :: half_grad_x(:,:), &
                        half_grad_y(:,:), &
                        half_h_grad_x(:,:), &
                        half_h_grad_y(:,:)
   integer :: i, j, order, sign 

   allocate(hx_tend(1-halo:num_cells+halo, 1-halo:num_cells+halo))
   allocate(hy_tend(1-halo:num_cells+halo, 1-halo:num_cells+halo))
   allocate(hx_tmp(1-halo:num_cells+halo, 1-halo:num_cells+halo))
   allocate(hy_tmp(1-halo:num_cells+halo, 1-halo:num_cells+halo))
   allocate(half_grad_x(0:num_cells, 1:num_cells))
   allocate(half_grad_y(1:num_cells, 0:num_cells))
   allocate(half_h_grad_x(0:num_cells, 1:num_cells))
   allocate(half_h_grad_y(1:num_cells, 0:num_cells))

    ! 1) Calculate even order laplacian
   hx_tmp(:,:) = h_old(:,:)
   do j = 1, num_cells
      do order = 1, diffusion_order / 2 -1
         do i = 0, num_cells+1
            hx_tend(i,j) = (hx_tmp(i+1,j) - 2 * hx_tmp(i,j) + hx_tmp(i-1,j)) / (dx**2)            
         end do
         if(order /= diffusion_order / 2 - 1) then
            hx_tmp(:,j) = hx_tend(:,j)
         end if
      end do
   end do

   hy_tmp(:,:) = h_old(:,:)
   do i = 1, num_cells
      do order = 1, diffusion_order / 2 -1
         do j = 0, num_cells+1
            hy_tend(i,j) = (hy_tmp(i,j+1) - 2 * hy_tmp(i,j) + hy_tmp(i,j-1)) / (dy**2)            
         end do
         if(order /= diffusion_order / 2 - 1) then
            hy_tmp(i,:) = hy_tend(i,:)
         end if
      end do
   end do 

   ! 2) Calculate gradient of laplcain above 
   do j = 1, num_cells
      do i = 0, num_cells
         half_grad_x(i,j) = (hx_tend(i+1,j) - hx_tend(i,j)) / dx
      end do 
   end do 
   do i = 1, num_cells
      do j = 0, num_cells
         half_grad_y(i,j) = (hy_tend(i,j+1) - hy_tend(i,j)) / dy
      end do 
   end do
   
   ! 3) Calculate gradient of the variable be diffused 
   do j = 1, num_cells
      do i = 0, num_cells
         half_h_grad_x(i,j) = (h_old(i+1,j) - h_old(i,j)) / dx
      end do 
   end do 
   do i = 1, num_cells
      do j = 0, num_cells
         half_h_grad_y(i,j) = (h_old(i,j+1) - h_old(i,j)) / dy
      end do 
   end do

   !4) decide the sign of diffusion flux to be zero or not
   sign = (-1)**(diffusion_order / 2 + 1)

   do j = 1, num_cells
      do i = 0, num_cells
         if (sign * half_grad_x(i,j) * half_h_grad_x(i,j) <= 0) then
            half_grad_x(i,j) = 0.0
         end if 
      end do 
   end do 
   do i = 1, num_cells
      do j = 0, num_cells
         if (sign * half_grad_y(i,j) * half_h_grad_y(i,j) <= 0) then
            half_grad_y(i,j) = 0.0
         end if
      end do 
   end do

   !5) Calculate the time tendency lastly.
   do j = 1, num_cells
      do i = 1, num_cells 
         hx_tend(i,j) = (half_grad_x(i,j) - half_grad_x(i-1,j)) / dx
      end do
   end do
   do i = 1, num_cells
      do j = 1, num_cells
         hy_tend(i,j) = (half_grad_y(i,j) - half_grad_y(i,j-1)) / dy
      end do 
   end do
   
   do j = 1, num_cells
      do i = 1, num_cells
        h_tend(i,j) = hx_tend(i,j) + hy_tend(i,j)
      end do 
   end do 
   
   return
   deallocate(hx_tend, hy_tend, hx_tmp, hy_tmp)
   deallocate(half_grad_x, half_grad_y, half_h_grad_x, half_h_grad_y)
end subroutine diffusion_limiter

subroutine full_periodic_boundary_condition(f, num_cells, halo)
   implicit none
   integer, intent(in) :: num_cells, halo
   real, intent(inout) :: f(1-halo:num_cells+halo, 1-halo:num_cells+halo)
   integer :: i, j
   do j = 1, num_cells
      do i = 1, halo
         f(1-i, j) = f(num_cells-i+1, j)
         f(num_cells+i, j) = f(1+i-1, j)
      end do 
   end do 
   do i = 1, num_cells
      do j = 1, halo
         f(i, 1-j) = f(i, num_cells-j+1)
         f(i, num_cells+j) = f(i, 1+j-1)
      end do
   end do 
end subroutine full_periodic_boundary_condition

end module diffusion_scheme


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