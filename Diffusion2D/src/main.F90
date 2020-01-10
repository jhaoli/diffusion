program main

   use namelist_mod
   use scheme_mod
   use history_mod

   implicit none

   real, allocatable :: h_old(:,:), h_new(:,:)
   real, allocatable :: x(:), y(:), h_tend(:,:)
   real :: diffusion_coef, dx,dy
   integer :: i, j, step,  order, sign, counter=0
   character(80) :: ofile
!!
   character(256) :: namelist_file_path
   
   call get_command_argument(1, namelist_file_path)
   if (namelist_file_path =='') then
     print*, 'You should give a namelist file path!'
     stop
   end if 

   call parse_namelist(namelist_file_path)
   

   ! Allocate spaces (we trade space/perforance for convenience)
   allocate(x(nx))
   allocate(y(ny))
   allocate(h_old(1-halo:nx+halo, 1-halo:ny+halo))
   allocate(h_new(1-halo:nx+halo, 1-halo:ny+halo))
  
   allocate(h_tend(1:nx, 1:ny))

   dx = (xmax - xmin) / nx
   dy = (ymax - ymin) / ny
   ! Set initial condition on uold, inititialize unew
   ! x = xmin + [((i-1)*dx,i=0,num_points+1)]
   ! y = ymin + [((j-1)*dy,j=0,num_points+1)]
   do i = 1, nx
      x(i) = (i-1)*dx
   end do 
   do j = 1, ny
      y(j) = (j-1)*dy
   end do 

   do j = 1, ny
     do i = 1, nx
       if (x(i)>= 0.3 .and. x(i) <=0.7 .and. y(j) >=0.3 .and. y(j)<=0.7) then
         h_old(i,j) = 0.4
       else
         h_old(i,j) = -0.4
       end if
     end do
   end do

   call full_periodic_boundary_condition(h_old)

   ! Output the initial values
   write(ofile,'(i3.3)') counter
   ofile = 'output'//trim(adjustl(ofile))//'.dat'
   open(11,file=ofile,status='unknown')
   print *, 'Ouput to ', ofile
   do j = 1, ny
      do i = 1, nx
      write(11,'(3f15.7)') x(i),y(j),h_old(i,j)
      end do
   end do
   close(11)

   call output(x,y,counter,h_old)

   stop
   counter = 1
   ! Perform time evolution process
   do step = 1, num_of_time_step
      if (diffusion_method == 'split') then
         call diffusion_split(h_old, h_tend, dx, dy)
      else if (diffusion_method == 'direct') then
         call diffusion_direct(h_old, h_tend, dx, dy)
      else if (diffusion_method == 'limiter') then
         call diffusion_limiter(h_old, h_tend, dx, dy)
      else
         print*, "no such method!"
      end if
 
      diffusion_coef = 1.0 / time_step_size * (dx/2.0)**diffusion_order
      sign = (-1)**(diffusion_order / 2 + 1)
      
      do j = 1, ny
         do i = 1, nx
         h_new(i,j) = h_old(i,j) + sign * time_step_size * diffusion_coef * h_tend(i,j)
         end do
      end do
      call full_periodic_boundary_condition(h_new)
      h_old(:,:) = h_new(:,:)
     
      ! Output the current solution
      if (mod(step, freq_output) == 0) then
         call output(x,y,counter, h_new)
         counter = counter + 1
      endif
   enddo

   ! Free spaces
   deallocate(x, y, h_old, h_new, h_tend)

end program main
