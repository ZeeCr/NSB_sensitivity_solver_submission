program incompressible_flow

  use param
  use aptofem_kernel
  use problem_options
  use fe_solution_restart_io
  use aptofem_read_control_file
  use aptofem_fe_solver_types
  use aptofem_linear_fe_solver
  use matrix_assembly_data_type

  use general_routines
  use geometry_routines
  use navier_stokes_jac_matrix_residual

  use coupled_data_type
  use fsi_aptofem_fe_solver
  use boundary_mapping

  use coupled_data_storage
  use boundary_data_storage
  use input_data_storage

  implicit none

  type(aptofem_keys), pointer :: aptofem_stored_keys

  character(len=aptofem_length_key_def) :: control_file, re_id
  character(len=aptofem_length_key_def) :: solve_type

  integer, dimension(:), allocatable :: bdry_map

  procedure(anal_soln_arg) :: anal_soln
  procedure(comp_bdry_no_arg) :: get_boundary_no

  procedure(d_bc_interpolate_arg) :: poiseuille_interpolate

  integer :: ierr

! Initialize package

  call AptoFEM_initialize(aptofem_stored_keys, 'aptofem_control_file.dat', './')

! Get user options

  call get_user_data('User_data', aptofem_stored_keys)

! Set up coupled scheme data

  call fsi_data_setup(scheme_coupled_data, problem_type='coupled', coupled_solns=coupled_solns)

! Set up bdry data

  call setup_bdry_storage()
  call setup_cvty_storage()
  call setup_otlt_storage()

! Set up input flow data, must be called after get_user_data

  call initialise_input_data(t_scaling)

! Create mesh

  ! Fluid mesh

  call create_mesh(fe_data_struct(1)%mesh_data, get_boundary_no, &
                   'Mesh_gen', aptofem_stored_keys)

  ! Set elements of mesh for cavities

  call set_element_cavity_no(fe_data_struct(1)%mesh_data)

! Get solver types

  call get_aptofem_key_definition('solver_type', &
                                  solve_type, trim(ns_solver_section), aptofem_stored_keys, ierr)

!Solver routines for Newton Navier-Stokes

  ! Fluid routines

  call store_subroutine_names(fe_data_struct(1)%fe_solver_subroutines, &
                              'assemble_jac_matrix_element', ns_jac_element_mat, 1)
  call store_subroutine_names(fe_data_struct(1)%fe_solver_subroutines, &
                              'assemble_jac_matrix_int_bdry_face', ns_dg_jac_face_mat, 1)
  call store_subroutine_names(fe_data_struct(1)%fe_solver_subroutines, &
                              'assemble_residual_element', ns_element_nonlin_residual, 1)
  call store_subroutine_names(fe_data_struct(1)%fe_solver_subroutines, &
                              'assemble_residual_int_bdry_face', ns_dg_face_nonlin_residual, 1)

! Create Solution

  ! Fluid solution

  call create_fe_solution(fe_data_struct(1)%soln_data, fe_data_struct(1)%mesh_data, trim(fe_solution_section), &
                          aptofem_stored_keys, anal_soln)

!Compute Solution

  allocate (bdry_map(no_bdries))
  call boundary_map(no_bdries, bdry_map)
  call coupled_data_setup(scheme_coupled_data, bdry_map)
  scheme_coupled_data%d_bc_interpolate_ptr => poiseuille_interpolate
  call scheme_coupled_data%d_bc_interpolate_ptr(scheme_coupled_data)

!Compute Solution

! BDF solver

  call fsi_compute_fe_solution_with_user_data(fe_data_struct, trim(ns_solver_section), aptofem_stored_keys, &
                                              scheme_coupled_data)

! Delete data structures

  deallocate (bdry_map)
  call delete_coupled_data(scheme_coupled_data)
  call delete_fe_data_structure(fe_data_struct, scheme_coupled_data)
  call delete_bdry_storage()
  call delete_cvty_storage()
  call delete_otlt_storage()

! Finalize package

  call AptoFEM_finalize(aptofem_stored_keys)

end program incompressible_flow
