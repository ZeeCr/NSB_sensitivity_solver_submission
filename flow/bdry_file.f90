!--------------------------------------------------------------------
subroutine get_boundary_no(boundary_no, strongly_enforced_bcs, global_point, &
                           face_coords, no_face_vert, problem_dim, mesh_data)
!--------------------------------------------------------------------
  use param
  use fe_mesh

  implicit none

  integer, intent(in) :: no_face_vert !< No of vertices of current face
  integer, intent(in) :: problem_dim !< Problem dimension
  real(db), dimension(no_face_vert, problem_dim), intent(in) :: face_coords
  !< Coordinates of the face
  real(db), dimension(problem_dim), intent(in) :: global_point
  !< Point on the boundary
  type(mesh), intent(in) :: mesh_data !< FE mesh
  integer, intent(inout) :: boundary_no !< Boundary identifier
  character(len=nvmax), intent(out) :: strongly_enforced_bcs !< Indicates
  !< which bcs are to be strongly imposed

! Local variables

  real(db) :: x, y, z, tol

  x = global_point(1)
  y = global_point(2)
  z = global_point(3)

  tol = 1.0d-7

  strongly_enforced_bcs = '0000'

end subroutine get_boundary_no

!--------------------------------------------------------------------
subroutine get_boundary_no_adv(boundary_no, strongly_enforced_bcs, global_point, &
                               face_coords, no_face_vert, problem_dim, mesh_data)
!--------------------------------------------------------------------
  use param
  use fe_mesh

  implicit none

  integer, intent(in) :: no_face_vert !< No of vertices of current face
  integer, intent(in) :: problem_dim !< Problem dimension
  real(db), dimension(no_face_vert, problem_dim), intent(in) :: face_coords
  !< Coordinates of the face
  real(db), dimension(problem_dim), intent(in) :: global_point
  !< Point on the boundary
  type(mesh), intent(in) :: mesh_data !< FE mesh
  integer, intent(inout) :: boundary_no !< Boundary identifier
  character(len=nvmax), intent(out) :: strongly_enforced_bcs !< Indicates
  !< which bcs are to be strongly imposed

  strongly_enforced_bcs = '0'

end subroutine get_boundary_no_adv

!--------------------------------------------------------------------
subroutine anal_soln(u, global_point, problem_dim, nv, bdry_no, t)
!--------------------------------------------------------------------
  use param
  use linear_algebra

  use boundary_data_storage
  use coupled_data_storage

  implicit none

  integer, intent(in) :: problem_dim !< problem dimension
  integer, intent(in) :: nv !< Number of variables in PDE system
  real(db), dimension(problem_dim), intent(in) :: global_point
  !< Point of evaluation.
  real(db), dimension(nv), intent(out) :: u !< Analytical solution
  integer, intent(in) :: bdry_no !< Boundary no
  real(db), intent(in) :: t !< Time
  real(db) :: x, y, z, tol

  real(db), dimension(problem_dim) :: n
  real(db), dimension(problem_dim) :: centre
  real(db) :: r
  real(db) :: radius

  integer :: bdry_index

  x = global_point(1)
  y = global_point(2)
  z = global_point(3)
  tol = 1.0d-9
  u = 0.0_db

  if (bdry_no > 1) then

    bdry_index = bdry_no_to_bdry_index(bdry_no)
    radius = bdry_storage(bdry_index)%radius
    centre = bdry_storage(bdry_index)%centre(1:problem_dim)
    n = bdry_storage(bdry_index)%normal(1:problem_dim)
    r = norm2(global_point - centre, problem_dim)

    u(1:problem_dim) = -scheme_coupled_data%peak_v(bdry_index)*n*abs(radius**2 - r**2)/radius**2

  end if

end subroutine anal_soln

!--------------------------------------------------------------------
subroutine neumann_bc(un, p, global_point, problem_dim, bdry_no)
!--------------------------------------------------------------------
  use param
  use problem_options
  use matrix_assembly_data_type
  implicit none

  integer, intent(in) :: problem_dim, bdry_no
  real(db), dimension(problem_dim), intent(out) :: un
  real(db), intent(out) :: p
  real(db), dimension(problem_dim), intent(in) :: global_point
  real(db) :: x, y, z, tol
  real(db), dimension(problem_dim + 1) :: u
  real(db), dimension(problem_dim + 1, problem_dim) :: grad_u
  real(db), dimension(problem_dim) :: normal
  tol = 1.0d-9

  x = global_point(1)
  y = global_point(2)
  z = global_point(3)

  un = 0.0_db

  p = 0.0_db
  !Traction (pressure) BC

end subroutine neumann_bc

!--------------------------------------------------------------------
subroutine calc_eff_permeability_reciprocal( &
  problem_dim, ele_no, point, perm, eff_perm)
!--------------------------------------------------------------------
  use param
  use linear_algebra

  use boundary_data_storage

  implicit none

  integer, intent(in) :: problem_dim, ele_no
  real(db), intent(in) :: perm !< Exterior permeability
  real(db), dimension(problem_dim) :: point
  real(db), intent(out) :: eff_perm

  integer :: outlet_no
  real(db) :: dist

  if (element_cavity_no(ele_no) == -1 .OR. element_outlet_no(ele_no) == -1) then
    dist = 0.0_db
  else if (element_cavity_no(ele_no) == -2) then
    outlet_no = element_outlet_no(ele_no)
    dist = norm2(point - otlt_storage(outlet_no)%interior_face_centre)
    ! CREATE A DOME
    if (dist > otlt_storage(outlet_no)%total_radius) then
      dist = dist - otlt_storage(outlet_no)%total_radius
    else
      dist = 0.0_db
    end if
  else if (element_outlet_no(ele_no) == -2) then
    ! In general_routines, set_element_cavity_no() already sets up the array
    ! so that element_cavity_no(ele_no) is the cavity no. of closest cavity
    ! to element ele_no
    !dist = cavity_ellipsoid_eval(problem_dim,element_cavity_no(ele_no),point)
    dist = &
      calc_dist_to_cavity_surface(problem_dim, element_cavity_no(ele_no), point)
  else
    call write_message(io_err, "ERROR: calc_eff_permeability_reciprocal")
    call write_message(io_err, "ERROR: case I didn't think of")
    write (io_err, *) element_cavity_no(ele_no), element_outlet_no(ele_no)
  end if

  eff_perm = perm*tanh(smoothing_freq*dist)

end subroutine calc_eff_permeability_reciprocal

!--------------------------------------------------------------------
subroutine forcing_fn(f, global_point, problem_dim, nv)
!--------------------------------------------------------------------
  use param
  implicit none

  integer, intent(in) :: problem_dim !< problem dimension
  integer, intent(in) :: nv !< Number of variables in PDE system
  real(db), dimension(nv), intent(out) :: f !< Forcing function
  real(db), dimension(problem_dim), intent(in) :: global_point
  !< Point of evaluation.

  f = 0.0_db

end subroutine forcing_fn
