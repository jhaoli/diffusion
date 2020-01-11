module mesh_mod
  
  use namelist_mod

  implicit none 

  private
  public mesh_type
  public mesh_init
  public mesh

  type mesh_type
    real, allocatable :: x(:), y(:)
    real :: dx, dy
  end type 
  
  type(mesh_type) mesh 

contains

  subroutine mesh_init()
    
    integer :: i, j

    allocate(mesh%x(nx+1))
    allocate(mesh%y(ny+1))

    mesh%dx = (xmax - xmin) / nx 
    mesh%dy = (ymax - ymin) / ny

    mesh%x = [( (i-1)*mesh%dx, i = 1, nx)]
    mesh%y = [( (j-1)*mesh%dy, j = 1, ny)]

  end subroutine mesh_init
end module mesh_mod