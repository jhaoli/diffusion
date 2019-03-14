
program diffusion
   implicit none
   real, allocatable :: h_old(:,:), h_new(:,:), hx_tmp(:,:), hy_tmp(:,:) 
   real, allocatable :: x(:), y(:), h_tend(:,:),hx_tend(:,:), hy_tend(:,:)
   real, allocatable :: half_grad_x(:,:), &
                        half_grad_y(:,:), &
                        half_h_grad_x(:,:), &
                        half_h_grad_y(:,:)
   real :: diffusion_coef, dx,dy, time_step_size, xmin, xmax, ymin, ymax
   integer :: num_points, num_of_time_step, freq_output, diffusion_order
   integer :: i, j, step, order, sign, counter=0
   character(80) :: ofile
   integer :: x_start, x_end, y_start, y_end

   namelist /params/ &
      time_step_size, &
      num_of_time_step, & 
      freq_output, & ! output frequency(every so many steps)
      xmin, &
      xmax, &
      ymin, &
      ymax, &
      num_points, &
      diffusion_order

   ! Input parameters
   open(10,file='namelist.input',status='old')
   read(10,nml=params)
   close(10)
  
   ! Allocate spaces (we trade space/perforance for convenience)
   allocate(x(0:num_points+1))
   allocate(y(0:num_points+1))
   allocate(h_old(0:num_points+1, 0:num_points+1))
   allocate(h_new(0:num_points+1, 0:num_points+1))
   allocate(hx_tmp(0:num_points+1, 0:num_points+1))
   allocate(hy_tmp(0:num_points+1, 0:num_points+1))
   allocate(h_tend(0:num_points+1, 0:num_points+1))
   allocate(hx_tend(0:num_points+1, 0:num_points+1))
   allocate(hy_tend(0:num_points+1, 0:num_points+1))

   allocate(half_grad_x(0:num_points+1, 0:num_points+1))
   allocate(half_grad_y(0:num_points+1, 0:num_points+1))
   allocate(half_h_grad_x(0:num_points+1, 0:num_points+1))
   allocate(half_h_grad_y(0:num_points+1, 0:num_points+1))
   
   dx = (xmax - xmin) / (num_points-1)
   dy = (ymax - ymin) / (num_points-1)
   ! Set initial condition on uold, inititialize unew
   x = xmin + [((i-1)*dx,i=0,num_points+1)]
   y = ymin + [((j-1)*dy,j=0,num_points+1)]

   h_old(:,:) = -0.4
   x_start = (0.28 - xmin) / dx + 1
   x_end = (0.72 - xmin) / dx + 1
   y_start = (0.28 - ymin) / dy + 1
   y_end = (0.72 - ymin) / dy + 1
   do j = y_start, y_end
      do i = x_start, x_end
         h_old(i,j) = 0.4
      end do
   end do
   ! Output the initial values
   write(ofile,'(i3.3)') counter
   ofile = 'output'//trim(adjustl(ofile))//'.dat'
   open(11,file=ofile,status='unknown')
   print *, 'Ouput to ', ofile
   do j = 1, num_points
      do i = 1, num_points
      write(11,'(3f15.7)') x(i),y(j),h_old(i,j)
      end do
   end do
   close(11)

   counter = 1
   

   ! Perform time evolution process
   do step = 1, num_of_time_step
      hx_tmp(:,:) = h_old(:,:)
      hy_tmp(:,:) = h_old(:,:)
         ! 1) Calculate even order laplacian
         do j = 0, num_points+1
            do order = 1, diffusion_order / 2 -1
               do i = 1, num_points
                  hx_tend(i,j) = (hx_tmp(i+1,j) - 2 * hx_tmp(i,j) + hx_tmp(i-1,j)) / (dx**2)            
               end do
               hx_tend(0, j) = 0
               hx_tend(num_points+1, j) = 0
               if(order /= diffusion_order / 2 - 1) then
                  hx_tmp(:,j) = hx_tend(:,j)
               end if
            end do
         end do

         do i = 0, num_points+1
            do order = 1, diffusion_order / 2 -1
               do j = 1, num_points
                  hy_tend(i,j) = (hy_tmp(i,j+1) - 2 * hy_tmp(i,j) + hy_tmp(i,j-1)) / (dy**2)            
               end do
               hy_tend(i, 0) = 0
               hy_tend(i, num_points+1) = 0
               if(order /= diffusion_order / 2 - 1) then
                  hy_tmp(i,:) = hy_tend(i,:)
               end if
            end do
         end do

         ! 2) Calculate gradient of 
         do j = 0, num_points + 1
            do i = 0, num_points
               half_grad_x(i,j) = (hx_tend(i+1,j) - hx_tend(i,j)) / dx
            end do 
         end do 
         do i = 0, num_points + 1
            do j = 0, num_points
               half_grad_y(i,j) = (hy_tend(i,j+1) - hy_tend(i,j)) / dy
            end do 
         end do
 
         ! 3) Calculate gradient of the variable be diffused 
         do j = 0, num_points + 1
            do i = 0, num_points
               half_h_grad_x(i,j) = (h_old(i+1,j) - h_old(i,j)) / dx
            end do 
         end do 
         do i = 0, num_points + 1
            do j = 0, num_points
               half_h_grad_y(i,j) = (h_old(i,j+1) - h_old(i,j)) / dy
            end do 
         end do

         !4) decide the sign of diffusion flux to be zero or not
         sign = (-1)**(diffusion_order / 2 + 1)
         do j = 0, num_points + 1
            do i = 0, num_points
               if (sign * half_grad_x(i,j) * half_h_grad_x(i,j) <=0) then
                  half_grad_x(i,j) = 0.0
               end if 
            end do 
         end do 
         do i = 0, num_points + 1
            do j = 0, num_points
               if (sign * half_grad_y(i,j) * half_h_grad_y(i,j) <=0) then
                  half_grad_y(i,j) = 0.0
               end if
            end do 
         end do

         !5) Calculate the time tendency finally.
         do j = 0, num_points+1
            do i = 1, num_points+1 
               hx_tend(i,j) = (half_grad_x(i,j) - half_grad_x(i-1,j)) / dx
            end do
         end do
         do i = 0, num_points + 1
            do j = 0, num_points
               hy_tend(i,j) = (half_grad_y(i,j) - half_grad_y(i,j-1)) / dy
            end do 
         end do

         h_tend(:,:) = hx_tend(:,:) + hy_tend(:,:)
               do j = 1, num_points
         do i = 1, num_points
              if( (h_tend(i,j) - 0.00) >= 0.001) then
               print*, i,j, h_tend(i,j)
            end if 
         end do 
      end do 
      stop
      diffusion_coef = 1.0 / time_step_size * (dx/2.0)**diffusion_order
      
      do j = 1, num_points
         do i = 1, num_points
         h_new(i,j) = h_old(i,j) + sign * time_step_size * diffusion_coef * h_tend(i,j)
         end do
      end do

      h_old(:,:) = h_new(:,:)
      ! Apply boundary conditions at two ends
      h_old(0,:) = -0.4      ! Left end point
      h_old(num_points+1,:) = -0.4      ! Right end point
      h_old(:,0) = -0.4 ! bottom end boundary
      h_old(:,num_points+1) = -0.4     ! up end boundary
      ! Output the current solution
      if (mod(step, freq_output) == 0) then
         write(ofile,'(i3.3)') counter
         ofile = 'output'//trim(adjustl(ofile))//'.dat'
         open(11,file=ofile,status='unknown')
         print *, 'Ouput to ', ofile
         do j = 1, num_points
            do i = 1, num_points
            write(11,'(3f15.7)') x(i),y(j), h_new(i,j)
            end do
         enddo
         close(11)
         counter = counter + 1
      endif
   enddo

   ! Free spaces
   deallocate(x, y, h_old,h_new,hx_tmp,hy_tmp, h_tend, hx_tend, hy_tend)
   deallocate(half_grad_x, half_grad_y, half_h_grad_x, half_h_grad_y)
end program diffusion
