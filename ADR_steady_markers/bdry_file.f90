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

  strongly_enforced_bcs = '0'

end subroutine get_boundary_no

!--------------------------------------------------------------------
subroutine anal_soln(u, global_point, problem_dim, nv, bdry_no, t)
!--------------------------------------------------------------------
  use param
  use linear_algebra
  use coupled_data_type

  implicit none

  integer, intent(in) :: problem_dim !< problem dimension
  integer, intent(in) :: nv !< Number of variables in PDE system
  real(db), dimension(problem_dim), intent(in) :: global_point
  !< Point of evaluation.
  real(db), dimension(nv), intent(out) :: u !< Analytical solution
  integer, intent(in) :: bdry_no !< Boundary no
  real(db), intent(in) :: t !< Time

  u = 0.0_db

end subroutine anal_soln

!--------------------------------------------------------------------
subroutine get_boundary_no_fluid(boundary_no, strongly_enforced_bcs, global_point, &
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

end subroutine get_boundary_no_fluid

!--------------------------------------------------------------------
subroutine anal_soln_fluid(u, global_point, problem_dim, nv, bdry_no, t)
!--------------------------------------------------------------------
  use param
  use linear_algebra
  use coupled_data_type

  implicit none

  integer, intent(in) :: problem_dim !< problem dimension
  integer, intent(in) :: nv !< Number of variables in PDE system
  real(db), dimension(problem_dim), intent(in) :: global_point
  !< Point of evaluation.
  real(db), dimension(nv), intent(out) :: u !< Analytical solution
  integer, intent(in) :: bdry_no !< Boundary no
  real(db), intent(in) :: t !< Time

  u = 0.0_db

end subroutine anal_soln_fluid
