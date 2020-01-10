
program diffusion
   implicit none
   real, allocatable :: w_old(:), w_new(:), w_tmp(:)
   real, allocatable :: x(:), w_tend(:)
   real :: diffusion_coef, dx, time_step_size, xmin, xmax
   integer :: num_cells, num_of_time_step, freq_output, diffusion_order
   integer :: i, step, halo, order, sign, counter=0
   character(80) :: ofile

   namelist /params/ &
      time_step_size, &
      num_of_time_step, & 
      freq_output, & ! output frequency(every so many steps)
      xmin, &
      xmax, &
      halo, &
      num_cells, &
      diffusion_order

   ! Input parameters
   open(10,file='namelist.input',status='old')
   read(10,nml=params)
   close(10)
  
   ! Allocate spaces (we trade space/perforance for convenience)
   allocate(x(num_cells+1))
   allocate(w_old(1-halo:num_cells+halo),w_new(1-halo:num_cells+halo), w_tmp(1-halo:num_cells+halo), w_tend(1-halo:num_cells+halo))
   dx = (xmax - xmin)/num_cells
   ! Set initial condition on uold, inititialize unew
   ! x = xmin + [((i-1)*dx,i=0,num_points+1)]

   ! where (0.28 <= x .and. x <= 0.72)
   !    w_old(:) = 0.4 
   ! else where
   !    w_old(:) = -0.4
   ! end where 
   do i = 1, num_cells+1
      x(i) = (i-1)*dx
   end do
   do i = 1, num_cells
      if(x(i)>=0.28 .and. x(i) <= 0.72) then
         w_old(i) = 0.4
      else
         w_old(i) = -0.4
      end if 
   end do 
   call full_boundary_condition(w_old)
   ! Output the initial values
   write(ofile,'(i3.3)') counter
   ofile = 'output'//trim(adjustl(ofile))//'.dat'
   open(11,file=ofile,status='unknown')
   print *, 'Ouput to ', ofile
   write(11,'(2f15.7)') (x(i),w_old(i),i=1,num_cells+1)
   close(11)

   counter = 1
   

   do step = 1, num_of_time_step
      w_tmp(:) = w_old(:)
      do order = 1, diffusion_order / 2
         do i = 1, num_cells
            w_tend(i) = (w_tmp(i+1) - 2 * w_tmp(i) + w_tmp(i-1)) / (dx**2)
         end do
         if(order /= diffusion_order / 2) then
            w_tmp(:) = w_tend(:)
         end if
      enddo
      

      diffusion_coef = 1.0 / time_step_size * (dx/2.0)**diffusion_order
      sign = (-1)**(diffusion_order / 2 + 1)
      
      do i = 1, num_cells
         w_new(i) = w_old(i) + sign * time_step_size * diffusion_coef * w_tend(i)
      end do

      
      ! Apply boundary conditions at two ends
      call full_periodic_boundary_condition(w_new)
      w_old(:) = w_new(:)
      ! Output the current solution
      if (mod(step, freq_output) == 0) then
         write(ofile,'(i3.3)') counter
         ofile = 'output'//trim(adjustl(ofile))//'.dat'
         open(11,file=ofile,status='unknown')
         print *, 'Ouput to ', ofile
         write(11,'(2f15.7)') (x(i),w_new(i),i=1,num_cells+1)
         close(11)
         counter = counter + 1
      endif
   enddo

   ! Free spaces
   deallocate(x,w_old,w_new,w_tmp,w_tend)

contains

subroutine full_periodic_boundary_condition(x)
   real, intent(inout) :: x(1-halo: num_cells+halo)
   integer :: i

   do i = 1, halo
      x(1-i) = x(num_cells-i+1)
      x(num_cells+i) = x(1+i-1)
   end do 
end subroutine full_periodic_boundary_condition

end program diffusion