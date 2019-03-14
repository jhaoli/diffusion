
program diffusion_flux_limiter
   implicit none
   real, allocatable :: w_init(:), w_old(:), w_new(:), w_tmp(:)
   real, allocatable :: x(:), w_tend(:)
   real, allocatable :: full_flux(:), half_grad(:), half_w_grad(:)
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
   allocate(full_flux(1-halo:num_cells+halo))
   allocate(half_grad(1-halo:num_cells+halo))
   allocate(half_w_grad(1-halo:num_cells+halo))

   dx = (xmax - xmin)/num_cells
   ! Set initial condition on uold, inititialize unew
   ! x = xmin + [((i-1)*dx,i=0,num_points+1)]
   ! where (0.24 <= x .and. x <= 0.72)
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
   call full_periodic_boundary_condition(w_old)
   ! Output the initial values
   write(ofile,'(i3.3)') counter
   ofile = 'output'//trim(adjustl(ofile))//'.dat'
   open(11,file=ofile,status='unknown')
   print *, 'Ouput to ', ofile
   write(11,'(2f15.7)') (x(i),w_old(i),i=1,num_cells+1)
   close(11)

   counter = 1
   
   ! Perform time evolution process
   do step = 1, num_of_time_step
      w_tmp = w_old
      do order = 1, diffusion_order / 2 -1
         ! 1) Calculate diffusion flux on cell center
         do i = 1, num_cells
            full_flux(i) = (w_tmp(i+1) - 2 * w_tmp(i) + w_tmp(i-1)) / (dx**2)
         end do
         call full_periodic_boundary_condition(full_flux)
         if(order /= diffusion_order / 2 -1) then
               w_tmp(:) = full_flux(:)
         end if
      end do

      ! 2) Calculate gradient of 2nd/4th order diffusion flux on cell boundary
      do i = 1, num_cells
         half_grad(i) = (full_flux(i+1) - full_flux(i)) / dx
      end do 
      
      
      ! 3) Calculate gradient of orign data on cell boundary
      do i = 1, num_cells
         half_w_grad(i) = (w_old(i+1) - w_old(i)) / dx
      end do
      ! half_w_grad(num_points+1) = 0

      ! 4) decide the sign on the cell bounday
      sign = (-1)**(diffusion_order / 2 + 1)
      do i = 1, num_cells
         if ( sign*half_grad(i) * half_w_grad(i) <= 0 ) then
            half_grad(i) = 0
         end if
      end do

      do i = 1, num_cells
         w_tend(i) = (half_grad(i) - half_grad(i-1)) / dx
      end do
     
     
!!!!      
      ! ! !1) cell boundary
      ! do i = 1, num_points-1
      !    full_flux(i) = (w_tmp(i+2) - 3 * w_tmp(i+1) + 3 * w_tmp(i) - w_tmp(i-1)) / (dx**3)
      ! end do
      ! full_flux(0) = 0
      ! full_flux(num_points) = 0
      ! full_flux(num_points+1) = 0

      ! do i = 0, num_points
      !     ! half_w_grad(i) =  (w_tmp(i+1) - w_tmp(i)) / dx
      !     half_w_grad(i) = ((w_tmp(i+1)+w_tmp(i))*0.5 - (w_tmp(i-1)+w_tmp(i))*0.5 ) / dx
      ! end do
      ! half_w_grad(num_points+1) = 0

      ! do i = 1, num_points-1
      !    if ( -1.0 * full_flux(i)*half_w_grad(i) <= 0 ) then
      !       full_flux(i) = 0
      !    end if
      ! end do

      ! do i = 0, num_points
      !    w_tend(i) = (full_flux(i) - full_flux(i-1))/ dx
      ! end do
   
      ! w_tend(num_points+1) = 0

      ! do i = 2, num_points-1
      !    w_tend(i) = (w_tmp(i+2) - 4 * w_tmp(i+1) + 6 * w_tmp(i) -4 * w_tmp(i-1) + w_tmp(i-2)) / (dx**4)
      ! end do


      diffusion_coef = 1.0 / time_step_size * (dx/2.0)**diffusion_order
      sign = (-1)**(diffusion_order / 2 + 1)
      
      do i = 1, num_cells
         w_new(i) = w_old(i) + sign * time_step_size * diffusion_coef * w_tend(i)
      end do
   
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
   deallocate(x, w_old, w_new,w_tmp,w_tend)
   deallocate(full_flux)
   deallocate(half_grad)
   deallocate(half_w_grad)
contains

subroutine full_periodic_boundary_condition(x)
   real, intent(inout) :: x(1-halo: num_cells+halo)
   integer :: i

   do i = 1, halo
      x(1-i) = x(num_cells-i+1)
      x(num_cells+i) = x(1+i-1)
   end do 
end subroutine full_periodic_boundary_condition

end program diffusion_flux_limiter
