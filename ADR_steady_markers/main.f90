program adr_markers

  use param
  use aptofem_kernel
  use problem_options
  use fe_solution_restart_io
  use aptofem_read_control_file
  use aptofem_fe_solver_types
  use aptofem_linear_fe_solver
  use matrix_assembly_data_type

  use coupled_data_type
  use fsi_aptofem_fe_solver

  use boundary_data_storage
  use coupled_data_storage
  use geometry_routines

  implicit none

  type(aptofem_keys), pointer :: aptofem_stored_keys
  type(coupled_user_data) :: scheme_coupled_data
  type(mesh) :: fluid_mesh
  type(solution) :: fluid_solution

  procedure(anal_soln_arg) :: anal_soln
  procedure(comp_bdry_no_arg) :: get_boundary_no
  procedure(anal_soln_arg) :: anal_soln_fluid
  procedure(comp_bdry_no_arg) :: get_boundary_no_fluid

  real(db), dimension(:), allocatable :: energy_markers

  integer :: i, j, ierr, no_energy_markers

  real(db) :: curr_t

  logical :: file_exists
  integer :: no_basal_veins, no_peripheral_veins, no_septal_veins, no_arteries, no_veins
  real(db), dimension(:), allocatable :: artery_radii, vein_radii, artery_areas, vein_areas
  ! Cotyledon nodal locations
  integer :: no_cells, max_nodes_per_cell
  integer, dimension(:), allocatable :: nodes_per_cell
  real(db), dimension(:, :, :), allocatable :: c_nodes

  integer :: id
  character(len=10) :: id_str
  character(len=aptofem_length_key_def) :: f_name_to_write

  character(len=aptofem_length_key_def) :: Re_str, Pe_str, DaRecipr_str

! Read integer from file
  open (unit=10, file='../data/id_idx.dat', status='old')
  read (10, *) id
  close (10)

  ! Get inlet / outlet info
  call read_no_outlets(no_basal_veins, no_peripheral_veins, no_septal_veins)
  call read_no_arteries(no_arteries)
  no_veins = no_basal_veins + no_peripheral_veins + no_septal_veins
  allocate (artery_radii(no_arteries), vein_radii(no_veins), &
            artery_areas(no_arteries), vein_areas(no_veins))
  artery_radii = 0.0_db
  vein_radii = 0.0_db
  artery_areas = 0.0_db
  vein_areas = 0.0_db
  call read_outlet_radii(no_basal_veins, no_peripheral_veins, no_septal_veins, &
                         vein_radii(1:no_basal_veins), &
                         vein_radii(1 + no_basal_veins:no_peripheral_veins + no_basal_veins), &
                        vein_radii(1 + no_peripheral_veins + no_basal_veins:no_septal_veins + no_peripheral_veins + no_basal_veins))
  call read_artery_radii(no_arteries, artery_radii)
  artery_areas(:) = pi*artery_radii(:)**2
  vein_areas(:) = pi*vein_radii(:)**2
  call read_cotyledon_node_locations(3, &
                                     no_cells, nodes_per_cell, max_nodes_per_cell, c_nodes)
  print *, "artery_areas:", artery_areas
  print *, "vein_areas:", vein_areas

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

! Create mesh

  ! Concentration mesh

  call create_mesh(fe_data_struct(1)%mesh_data, get_boundary_no, &
                   'Mesh_gen', aptofem_stored_keys)

  ! Fluid mesh

  call create_mesh(fluid_mesh, get_boundary_no_fluid, &
                   'Mesh_gen', aptofem_stored_keys)

! Set elements of mesh for cavities

  call set_element_cavity_no(fe_data_struct(1)%mesh_data)

!Solver routines for Newton Navier-Stokes

! Create Solution

  ! Concentration solution

  call create_fe_solution(fe_data_struct(1)%soln_data, fe_data_struct(1)%mesh_data, 'FE_Adv_Solution', &
                          aptofem_stored_keys, anal_soln)

  ! Fluid solution

  call create_fe_solution(fluid_solution, fluid_mesh, 'FE_Solution_full_system', &
                          aptofem_stored_keys, anal_soln_fluid)

! Read Solution and calculate energy_markers
  no_energy_markers = no_veins + 2*(no_cells + 1) + 5
  allocate (energy_markers(no_energy_markers))

  write (Re_str, '(es16.6)') Re
  write (Pe_str, '(es16.4)') Pe
  write (DaRecipr_str, '(es16.4)') Da_reciprocal
  write (f_name_to_write, '(i0, a)') &
    id, '_Re'//trim(adjustl(Re_str))// &
    '_Pe'//trim(adjustl(Pe_str))// &
    '_DaRecipr'//trim(adjustl(DaRecipr_str))//'_adr_steady_markers.csv'

! Save markers info into CSV format, create file if doesn't exist, else only open
  inquire (file=trim(adjustl(f_name_to_write)), exist=file_exists)
  if (.NOT. file_exists) then
    open (unit=104, file=trim(adjustl(f_name_to_write)), status='new')
    write (104, '("Dm")', advance='no')
    do j = 1, no_basal_veins
      write (104, '(a2,i0)', advance='no') ",b", j
    end do
    do j = 1, no_peripheral_veins
      write (104, '(a2,i0)', advance='no') ",p", j
    end do
    do j = 1, no_septal_veins
      write (104, '(a2,i0)', advance='no') ",s", j
    end do
    do j = 1, no_cells
      write (104, '(a7,i0)', advance='no') ",uptake", j
    end do
    write (104, '(a13)', advance='no') ",total_uptake"
    do j = 1, no_cells
      write (104, '(a14,i0)', advance='no') ",scaled_uptake", j
    end do
    write (104, '(a)', advance='no') ",scaled_total_uptake"
    write (104, '(a)', advance='no') ",adv_flux_in"
    write (104, '(a)', advance='no') ",adv_flux_out"
    write (104, '(a)', advance='no') ",total_surface_flux"
    write (104, '(a)', advance='no') ",in_out_total_flux"
    write (104, '(a)', advance='no') ",in_out_flux_bc_adhere"
    write (104, *)
  else
    open (unit=104, file=trim(adjustl(f_name_to_write)), &
          status='old', position="append", action="write")
  end if

  ! Concentration solution

  call read_solution_for_restart(fe_data_struct(1)%mesh_data, fe_data_struct(1)%soln_data, 2, &
                                 'adv_steady', 2, '../ADR_steady/output/restart/')

  ! Fluid solution

  call read_solution_for_restart(fluid_mesh, fluid_solution, &
                                 81, 'soln_1', 2, '../flow/output/restart/')

  ! Markers

  call energy_marker(fe_data_struct(1)%mesh_data, fe_data_struct(1)%soln_data, &
                     no_arteries, no_veins, &
                     fluid_mesh, fluid_solution, &
                     no_cells, nodes_per_cell, max_nodes_per_cell, c_nodes, &
                     no_energy_markers, energy_markers)

  write (104, '(f0.8)', advance='no') Dm
  do j = 1, no_energy_markers
    write (104, '(a1,f0.9)', advance='no') ',', energy_markers(j)
  end do
  write (104, *)

  close (unit=104)
  deallocate (energy_markers)

! Delete data structures

  call delete_coupled_data(scheme_coupled_data)
  call delete_fe_data_structure(fe_data_struct, scheme_coupled_data)
  call delete_bdry_storage()
  call delete_cvty_storage()
  call delete_otlt_storage()

  deallocate (artery_radii, vein_radii, artery_areas, vein_areas, nodes_per_cell, c_nodes)

! Finalize package

  call AptoFEM_finalize(aptofem_stored_keys)

end program adr_markers

!--------------------------------------------------------------------
subroutine energy_marker(mesh_data, soln_data, &
                         no_arteries, no_veins, &
                         fluid_mesh, fluid_soln, &
                         no_cells, nodes_per_cell, max_nodes_per_cell, c_nodes, &
                         no_energy_markers, energy_markers)
!--------------------------------------------------------------------
  use param
  use fe_mesh
  use fe_solution
  use basis_fns_storage_type
  use aptofem_fe_matrix_assembly
  use problem_options

  use geometry_routines

  implicit none

  type(mesh), intent(inout) :: mesh_data, fluid_mesh !*FD FE mesh
  type(solution), intent(inout) :: soln_data, fluid_soln !*FD FE solution
  integer, intent(in) :: no_arteries, no_veins, no_energy_markers
  real(db), dimension(no_energy_markers), intent(out) :: energy_markers

  integer, intent(in) :: no_cells, max_nodes_per_cell
  integer, dimension(no_cells), intent(in) :: nodes_per_cell
  real(db), dimension(mesh_data%problem_dim, max_nodes_per_cell, no_cells), &
    intent(in) :: c_nodes

! Local variables

  type(basis_storage) :: fe_basis_info
  character(len=aptofem_length_key_def) :: control_parameter
  integer :: no_eles, no_nodes, no_faces, problem_dim, no_pdes, &
             i, j, k, qk, iv, no_quad_points, npinc, &
             no_quad_points_volume_max, no_quad_points_face_max, &
             bdry_face, dim_soln_coeff, no_pdes_fluid
  real(db), dimension(:, :), allocatable :: global_points_ele
  real(db), dimension(:), allocatable :: quad_weights_ele, jacobian
  real(db), dimension(:, :), allocatable :: gradient_u
  real(db), dimension(:), allocatable :: u, adv_vel
  real(db) :: l2_norm, h1_semi_norm, dg_norm, full_dispenal, h2_norm, hessian_norm
  real(db), dimension(:), allocatable :: quad_weights_face, face_jacobian, &
                                         dispenal
  real(db), dimension(:, :), allocatable :: global_points_face, face_normals, hessian_uh
  integer, dimension(:, :), allocatable :: global_dof_numbers1, global_dof_numbers2, &
                                           global_dof_numbers
  integer, dimension(2) :: loc_face_no, neighbors
  integer, dimension(:), allocatable :: no_dofs_per_variable1, &
                                        no_dofs_per_variable2, no_dofs_per_variable
  real(db), dimension(:, :, :), allocatable :: hessian_u

  integer, dimension(max_nodes_per_cell) :: sgn_vec ! Indicates whether node is on particular side of edge
  real(db) :: flux_in, flux_out, vol, eff_Dm, tol, pt_eval, &
              total_uptake, scaled_total_uptake, total_surface_flux, &
              in_out_total_flux, in_out_flux_bc_adhere
  real(db), dimension(mesh_data%problem_dim - 1) :: edge_dir, pt_vec
  real(db), dimension(no_cells) :: cotyledon_vol, cotyledon_uptake, scaled_cotyledon_uptake
  real(db), dimension(no_energy_markers) :: energy_marker_in, energy_marker_out

  real(db), dimension(no_arteries) :: face_flux_in
  real(db), dimension(no_veins) :: face_area_out, face_flux_out

  dim_soln_coeff = get_dim_soln_coeff(soln_data)
  no_pdes = get_no_pdes(soln_data)

  no_pdes_fluid = get_no_pdes(fluid_soln)

  call get_mesh_info(no_eles, no_nodes, no_faces, problem_dim, &
                     mesh_data)

  npinc = 4
  call compute_max_no_quad_points(no_quad_points_volume_max, &
                                  no_quad_points_face_max, mesh_data, soln_data, npinc)

  control_parameter = 'fo_deriv_uh_ele'
  call initialize_fe_basis_storage(fe_basis_info, control_parameter, soln_data, &
                                   problem_dim, no_quad_points_volume_max, no_quad_points_face_max)

  allocate (gradient_u(no_pdes, problem_dim))
  allocate (u(no_pdes))
  allocate (hessian_uh(problem_dim, problem_dim))
  allocate (hessian_u(no_pdes, problem_dim, problem_dim))
  allocate (no_dofs_per_variable(dim_soln_coeff))
  allocate (global_points_ele(problem_dim, no_quad_points_volume_max))
  allocate (jacobian(no_quad_points_volume_max))
  allocate (quad_weights_ele(no_quad_points_volume_max))
  allocate (global_dof_numbers(dim_soln_coeff, no_ele_dofs_per_var_max))
  allocate (adv_vel(no_pdes_fluid))

  face_area_out = 0.0_db
  face_flux_in = 0.0_db
  face_flux_out = 0.0_db
  flux_in = 0.0_db
  flux_out = 0.0_db
  total_surface_flux = 0.0_db
  in_out_total_flux = 0.0_db
  in_out_flux_bc_adhere = 0.0_db

  vol = 0.0_db
  total_uptake = 0.0_db
  cotyledon_vol = 0.0_db
  cotyledon_uptake = 0.0_db

  tol = 1.0d-6

  do k = 1, no_eles

    call element_integration_info(dim_soln_coeff, problem_dim, mesh_data, &
                                  soln_data, k, npinc, no_quad_points_volume_max, &
                                  no_quad_points, global_points_ele, jacobian, quad_weights_ele, &
                                  global_dof_numbers, no_dofs_per_variable, fe_basis_info)

    do qk = 1, no_quad_points

      call eval_tanh_smooth( &
        problem_dim, k, global_points_ele(:, qk), Dm, eff_Dm)

      u = uh_element(fe_basis_info, no_pdes, qk)

      total_uptake = total_uptake + jacobian(qk)*quad_weights_ele(qk) &
                     *eff_Dm*u(1)

      ! Only consider elements with uptake
      if (eff_Dm > tol) then
        vol = vol + jacobian(qk)*quad_weights_ele(qk)

        do i = 1, no_cells
          sgn_vec = 0
          do j = 1, nodes_per_cell(i)
            if (j /= nodes_per_cell(i)) then
              edge_dir = c_nodes(1:problem_dim - 1, j + 1, i) - c_nodes(1:problem_dim - 1, j, i)
            else
              edge_dir = c_nodes(1:problem_dim - 1, 1, i) - c_nodes(1:problem_dim - 1, j, i)
            end if
            pt_vec = global_points_ele(1:problem_dim - 1, qk) - c_nodes(1:problem_dim - 1, j, i)

            pt_eval = edge_dir(1)*pt_vec(2) - edge_dir(2)*pt_vec(1)

            if (pt_eval > tol) then
              sgn_vec(j) = 1
            else
              sgn_vec(j) = -1
            end if

            if (j > 1) then
              if (sgn_vec(j) /= sgn_vec(j - 1)) then
                exit
              end if
            end if

          end do

          if (SUM(sgn_vec) == nodes_per_cell(i)) then
            cotyledon_vol(i) = cotyledon_vol(i) + jacobian(qk)*quad_weights_ele(qk)
            cotyledon_uptake(i) = cotyledon_uptake(i) + jacobian(qk)*quad_weights_ele(qk)*eff_Dm*u(1)
            exit
          end if

        end do
      end if

    end do

  end do

  scaled_total_uptake = total_uptake/vol
  scaled_cotyledon_uptake = cotyledon_uptake/cotyledon_vol

  call delete_fe_basis_storage(fe_basis_info)

  deallocate (no_dofs_per_variable, global_points_ele, jacobian, &
              quad_weights_ele, global_dof_numbers, hessian_u, hessian_uh)

  control_parameter = 'fo_deriv_uh_face'
  call initialize_fe_basis_storage(fe_basis_info, control_parameter, soln_data, &
                                   problem_dim, no_quad_points_volume_max, no_quad_points_face_max)

  allocate (global_points_face(problem_dim, no_quad_points_face_max))
  allocate (face_jacobian(no_quad_points_face_max))
  allocate (face_normals(problem_dim, no_quad_points_face_max))
  allocate (quad_weights_face(no_quad_points_face_max))
  allocate (no_dofs_per_variable1(dim_soln_coeff))
  allocate (no_dofs_per_variable2(dim_soln_coeff))
  allocate (dispenal(dim_soln_coeff))
  allocate (global_dof_numbers1(dim_soln_coeff, no_ele_dofs_per_var_max))
  allocate (global_dof_numbers2(dim_soln_coeff, no_ele_dofs_per_var_max))

  do k = 1, no_faces

    call face_integration_info(dim_soln_coeff, problem_dim, mesh_data, soln_data, &
                               k, neighbors, loc_face_no, npinc, no_quad_points_face_max, &
                               no_quad_points, global_points_face, face_jacobian, face_normals, &
                               quad_weights_face, global_dof_numbers1, no_dofs_per_variable1, &
                               bdry_face, global_dof_numbers2, no_dofs_per_variable2, &
                               fe_basis_info)

! Boundary face

    do qk = 1, no_quad_points

      call compute_uh_glob_pt(adv_vel(:), &
                              no_pdes_fluid, neighbors(1), &
                              global_points_face(:, qk), problem_dim, fluid_mesh, fluid_soln)

      u = uh_face1(fe_basis_info, no_pdes, qk)
      gradient_u(no_pdes, :) = grad_uh_face1(fe_basis_info, problem_dim, no_pdes, &
                                             qk, 1)

      if (bdry_face > 0) then

        total_surface_flux = total_surface_flux + &
                             face_jacobian(qk)*quad_weights_face(qk)*dot_product( &
                             (1.0/Pe)*gradient_u(no_pdes, 1:problem_dim) - u(1)*adv_vel(1:problem_dim), &
                             face_normals(:, qk))

      end if

      if (bdry_face > 1) then

        in_out_total_flux = in_out_total_flux + &
                            face_jacobian(qk)*quad_weights_face(qk)*dot_product( &
                            (1.0/Pe)*gradient_u(no_pdes, 1:problem_dim) - u(1)*adv_vel(1:problem_dim), &
                            face_normals(:, qk))

      end if

      if (bdry_face > 1 .AND. bdry_face < 1 + no_arteries + 1) then

        face_flux_in(bdry_face - 1) = face_flux_in(bdry_face - 1) &
                                      - face_jacobian(qk)*quad_weights_face(qk)*u(1)* &
                                      dot_product(adv_vel(1:problem_dim), face_normals(:, qk))

        in_out_flux_bc_adhere = in_out_flux_bc_adhere + &
                                face_jacobian(qk)*quad_weights_face(qk)*dot_product( &
                                (1.0/Pe)*gradient_u(no_pdes, 1:problem_dim) - adv_vel(1:problem_dim), &
                                face_normals(:, qk))

      end if

      if (bdry_face > 100 .AND. bdry_face < 100 + no_veins + 1) then

        face_area_out(bdry_face - 100) = face_area_out(bdry_face - 100) &
                                         + face_jacobian(qk)*quad_weights_face(qk)

        face_flux_out(bdry_face - 100) = face_flux_out(bdry_face - 100) &
                                         + face_jacobian(qk)*quad_weights_face(qk)*u(1)* &
                                         dot_product(adv_vel(1:problem_dim), face_normals(:, qk))

        in_out_flux_bc_adhere = in_out_flux_bc_adhere + &
                                face_jacobian(qk)*quad_weights_face(qk)*dot_product( &
                                -u(1)*adv_vel(1:problem_dim), &
                                face_normals(:, qk))

      end if

    end do

  end do

  deallocate (gradient_u, global_points_face, face_jacobian, face_normals, quad_weights_face, &
              no_dofs_per_variable1, no_dofs_per_variable2, dispenal, &
              global_dof_numbers1, global_dof_numbers2, u, adv_vel)

  call delete_fe_basis_storage(fe_basis_info)

  flux_in = SUM(face_flux_in)
  flux_out = SUM(face_flux_out)

  energy_markers(1:no_veins) = face_flux_out(:)/flux_out
  do i = 1, no_cells
    energy_markers(no_veins + i) = cotyledon_uptake(i)
    energy_markers(no_veins + no_cells + 1 + i) = scaled_cotyledon_uptake(i)
  end do
  energy_markers(no_veins + no_cells + 1) = total_uptake
  energy_markers(no_veins + 2*(no_cells + 1)) = scaled_total_uptake
  energy_markers(no_veins + 2*(no_cells + 1) + 1) = flux_in
  energy_markers(no_veins + 2*(no_cells + 1) + 2) = flux_out
  energy_markers(no_veins + 2*(no_cells + 1) + 3) = total_surface_flux
  energy_markers(no_veins + 2*(no_cells + 1) + 4) = in_out_total_flux
  energy_markers(no_veins + 2*(no_cells + 1) + 5) = in_out_flux_bc_adhere

end subroutine energy_marker
