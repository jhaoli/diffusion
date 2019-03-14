
program diffusion_main
   use diffusion_scheme

   implicit none
   real, allocatable :: h_old(:,:), h_new(:,:)
   real, allocatable :: x(:), y(:), h_tend(:,:)
   real :: diffusion_coef, dx,dy, time_step_size, xmin, xmax, ymin, ymax
   integer :: num_cells, num_of_time_step, freq_output, diffusion_order
   integer :: i, j, step, halo, order, sign, counter=0
   character(80) :: ofile, diffusion_method
   integer :: x_start, x_end, y_start, y_end

   namelist /params/ &
      time_step_size, &
      num_of_time_step, & 
      freq_output, & ! output frequency(every so many steps)
      xmin, &
      xmax, &
      ymin, &
      ymax, &
      halo, &
      num_cells, &
      diffusion_order, &
      diffusion_method

   ! Input parameters
   open(10,file='namelist.input',status='old')
   read(10,nml=params)
   close(10)
  
   ! Allocate spaces (we trade space/perforance for convenience)
   allocate(x(num_cells+1))
   allocate(y(num_cells+1))
   allocate(h_old(1-halo:num_cells+halo, 1-halo:num_cells+halo))
   allocate(h_new(1-halo:num_cells+halo, 1-halo:num_cells+halo))
  
   allocate(h_tend(1:num_cells, 1:num_cells))

   dx = (xmax - xmin) / num_cells
   dy = (ymax - ymin) / num_cells
   ! Set initial condition on uold, inititialize unew
   ! x = xmin + [((i-1)*dx,i=0,num_points+1)]
   ! y = ymin + [((j-1)*dy,j=0,num_points+1)]
   do i = 1, num_cells
      x(i) = (i-1)*dx
   end do 
   do j = 1, num_cells
      y(j) = (j-1)*dy
   end do 
   h_old(:,:) = -0.4
   x_start = ((xmax-xmin)/3 - xmin) / dx + 1
   x_end = ((xmax-xmin)*2/3 - xmin) / dx + 1
   y_start = ((ymax-ymin)/3 - ymin) / dy + 1
   y_end = ((ymax-ymin)*2/3 - ymin) / dy + 1
   do j = y_start, y_end
      do i = x_start, x_end
         h_old(i,j) = 0.4
      end do
   end do
   call full_periodic_boundary_condition(h_old, num_cells, halo)
   ! Output the initial values
   write(ofile,'(i3.3)') counter
   ofile = 'output'//trim(adjustl(ofile))//'.dat'
   open(11,file=ofile,status='unknown')
   print *, 'Ouput to ', ofile
   do j = 1, num_cells
      do i = 1, num_cells
      write(11,'(3f15.7)') x(i),y(j),h_old(i,j)
      end do
   end do
   close(11)

   counter = 1
   

   ! Perform time evolution process
   do step = 1, num_of_time_step
      if (diffusion_method == 'split') then
         call diffusion_split(h_old, h_tend, num_cells, diffusion_order, dx, dy, halo)
      else if (diffusion_method == 'direct') then
         call diffusion_direct(h_old, h_tend, num_cells, diffusion_order, dx, dy, halo)
      else if (diffusion_method == 'limiter') then
         call diffusion_limiter(h_old, h_tend, num_cells, diffusion_order, dx, dy, halo)
      else
         print*, "no such method!"
      end if
 
      diffusion_coef = 1.0 / time_step_size * (dx/2.0)**diffusion_order
      sign = (-1)**(diffusion_order / 2 + 1)
      
      do j = 1, num_cells
         do i = 1, num_cells
         h_new(i,j) = h_old(i,j) + sign * time_step_size * diffusion_coef * h_tend(i,j)
         end do
      end do
      call full_periodic_boundary_condition(h_new, num_cells, halo)
      h_old(:,:) = h_new(:,:)
     
      ! Output the current solution
      if (mod(step, freq_output) == 0) then
         write(ofile,'(i3.3)') counter
         ofile = 'output'//trim(adjustl(ofile))//'.dat'
         open(11,file=ofile,status='unknown')
         print *, 'Ouput to ', ofile
         do j = 1, num_cells
            do i = 1, num_cells
            write(11,'(3f15.7)') x(i),y(j), h_new(i,j)
            end do
         enddo
         close(11)
         counter = counter + 1
      endif
   enddo

   ! Free spaces
   deallocate(x, y, h_old, h_new, h_tend)

end program diffusion_main
