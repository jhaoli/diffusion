module mesh_mod

  use namelist_mod

  implicit none 
  
  private 

  public mesh_init
  public mesh_final
  public mesh_type
  public mesh 

  type mesh_type
    real, allocatable :: x(:)
    real dx  
  end type mesh_type

  type(mesh_type)  mesh

contains 

  subroutine mesh_init()
    
    integer :: i

    allocate(mesh%x(1:nx))

    mesh%dx = (xmax - xmin)/nx
    
    mesh%x = [(xmin + mesh%dx * 0.5 + (i-1)*mesh%dx, i = 1, nx)]

  end subroutine mesh_init

  subroutine mesh_final()
    ! class(mesh_type), intent(inout) :: this

    if(allocated(mesh%x)) deallocate(mesh%x)

  end subroutine mesh_final

end module mesh_mod