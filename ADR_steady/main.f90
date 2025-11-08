program ADR_markers

  use param
  use aptofem_kernel
  use problem_options
  use fe_solution_restart_io
  use aptofem_read_control_file
  use aptofem_fe_solver_types
  use aptofem_linear_fe_solver
  use matrix_assembly_data_type

  use adr_av_jac_matrix_residual
  use geometry_routines

  use coupled_data_type

  use coupled_data_storage
  use boundary_data_storage

  implicit none

  type(aptofem_keys), pointer :: aptofem_stored_keys
  type(solution) :: dg_soln_data

  character(len=aptofem_length_key_def) :: control_file
  character(len=aptofem_length_key_def) :: solve_type

  procedure(anal_soln_arg) :: anal_soln
  procedure(anal_soln_arg) :: anal_soln_fluid
  procedure(comp_bdry_no_arg) :: get_boundary_no
  procedure(comp_bdry_no_arg) :: get_boundary_no_fluid

  integer :: ierr
  real(db) :: Pe_final

! Initialize package

  control_file = 'aptofem_control_file.dat'
  call AptoFEM_initialize(aptofem_stored_keys, trim(adjustl(control_file)), './')

! Get user options and store Pe for continuation

  call get_user_data('User_data', aptofem_stored_keys)
  Pe_final = Pe

! Set up coupled scheme data

  call fsi_data_setup(scheme_coupled_data, problem_type='coupled', coupled_solns=coupled_solns)

! Set up bdry data

  call setup_bdry_storage()
  call setup_cvty_storage()
  call setup_otlt_storage()

! Create mesh

  ! Advection mesh

  call create_mesh(fe_data_struct(1)%mesh_data, get_boundary_no, &
                   'Mesh_gen', aptofem_stored_keys)

  ! External mesh

  call create_mesh(fe_data_struct(1)%external_mesh, get_boundary_no_fluid, &
                   'Mesh_gen', aptofem_stored_keys)

  ! Set elements of mesh for cavities

  call set_element_cavity_no(fe_data_struct(1)%mesh_data)

!Solver routines for Newton Navier-Stokes

  ! Advection routines

  call store_subroutine_names(fe_data_struct(1)%fe_solver_subroutines, &
                              'assemble_jac_matrix_element', adr_av_jac_element_mat, 1)
  call store_subroutine_names(fe_data_struct(1)%fe_solver_subroutines, &
                              'assemble_jac_matrix_int_bdry_face', adr_av_dg_jac_face_mat, 1)

  call store_subroutine_names(fe_data_struct(1)%fe_solver_subroutines, &
                              'assemble_residual_element', adr_av_element_nonlin_residual, 1)
  call store_subroutine_names(fe_data_struct(1)%fe_solver_subroutines, &
                              'assemble_residual_int_bdry_face', adr_av_dg_face_nonlin_residual, 1)

! Create Solution

  ! Advection solution

  call create_fe_solution(fe_data_struct(1)%soln_data, fe_data_struct(1)%mesh_data, trim(adv_solution_section), &
                          aptofem_stored_keys, anal_soln)

  ! DG Fluid solution

  call create_fe_solution(fe_data_struct(1)%external_soln, fe_data_struct(1)%external_mesh, trim(fe_solution_section), &
                          aptofem_stored_keys, anal_soln_fluid)

  call read_solution_for_restart(fe_data_struct(1)%external_mesh, fe_data_struct(1)%external_soln, &
                                 81, 'soln_1', 2, '../flow/output/restart/')

  ! Continuation
  Pe = Pe_final/1000.0
  ierr = replace_key_val_acf_real_db(control_file, 'newton_tolerance', 1.0d-5)

  call compute_fe_solution(fe_data_struct(1)%soln_data, &
                           fe_data_struct(1)%mesh_data, fe_data_struct(1)%fe_solver_subroutines, &
                           trim(solver_section), aptofem_stored_keys, scheme_coupled_data)

  ! Final value
  Pe = Pe_final
  ierr = replace_key_val_acf_real_db(control_file, 'newton_tolerance', 1.0d-8)

  call compute_fe_solution(fe_data_struct(1)%soln_data, &
                           fe_data_struct(1)%mesh_data, fe_data_struct(1)%fe_solver_subroutines, &
                           trim(solver_section), aptofem_stored_keys, scheme_coupled_data)

  call write_fe_data('Output_mesh_solution', aptofem_stored_keys, &
                     1, fe_data_struct(1)%mesh_data, fe_data_struct(1)%soln_data)

  call write_solution_for_restart( &
    fe_data_struct(1)%soln_data, fe_data_struct(1)%mesh_data, aptofem_run_number, &
    'adv_steady', './output/restart/')

! Delete data structures

  !call delete_solution(dg_soln_data)
  call delete_coupled_data(scheme_coupled_data)
  call delete_fe_data_structure(fe_data_struct, scheme_coupled_data)
  call delete_bdry_storage()
  call delete_cvty_storage()
  call delete_otlt_storage()

! Finalize package

  call AptoFEM_finalize(aptofem_stored_keys)

end program ADR_markers
