program flow_markers

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

  procedure(anal_soln_arg) :: anal_soln
  procedure(comp_bdry_no_arg) :: get_boundary_no
  
  integer :: i,j,ierr
  integer, parameter :: no_files = 71
  integer, dimension(no_files) :: file_no_arr
  real(db) :: curr_t
  
  integer, parameter :: no_energy_markers = 58
  real(db), dimension(no_energy_markers,no_files+1) :: energy_markers

  integer :: no_basal_veins,no_peripheral_veins,no_septal_veins,no_arteries,no_veins
  real(db), dimension(:), allocatable :: artery_radii,vein_radii,artery_areas,vein_areas

  integer :: id
  character(len=10) :: id_str
  character(len=aptofem_length_key_def) :: f_name_to_write

  character(len=aptofem_length_key_def) :: Re_str,DaRecipr_str

! Read integer from file
  open(unit=10, file='../data/id_idx.dat', status='old')
  read(10, *) id
  close(10)

  energy_markers = 0.0_db
  do i = 1,no_files
	  file_no_arr(i) = 70+(i-1)*1
  end do

  ! Get inlet / outlet info
  call read_no_outlets(no_basal_veins,no_peripheral_veins,no_septal_veins)
  call read_no_arteries(no_arteries)
  no_veins = no_basal_veins+no_peripheral_veins+no_septal_veins
  allocate(artery_radii(no_arteries),vein_radii(no_veins), &
    artery_areas(no_arteries),vein_areas(no_veins))
  artery_radii = 0.0_db
  vein_radii = 0.0_db
  artery_areas = 0.0_db
  vein_areas = 0.0_db
  call read_outlet_radii(no_basal_veins,no_peripheral_veins,no_septal_veins, &
    vein_radii(1:no_basal_veins), &
    vein_radii(1+no_basal_veins:no_peripheral_veins+no_basal_veins), &
    vein_radii(1+no_peripheral_veins+no_basal_veins:no_septal_veins+no_peripheral_veins+no_basal_veins))
  call read_artery_radii(no_arteries,artery_radii)
  artery_areas(:) = pi*artery_radii(:)**2
  vein_areas(:) = pi*vein_radii(:)**2

! Initialize package

  call AptoFEM_initialize(aptofem_stored_keys,'aptofem_control_file.dat','./')

! Get user options

  call get_user_data('User_data',aptofem_stored_keys)
  
! Set up coupled scheme data

  call fsi_data_setup(scheme_coupled_data,problem_type='coupled',coupled_solns=coupled_solns)
  
! Set up bdry data
  
  call setup_bdry_storage()
  call setup_cvty_storage()
  call setup_otlt_storage()
  
! Create mesh

  ! Fluid mesh

  call create_mesh(fe_data_struct(1)%mesh_data,get_boundary_no, &
		'Mesh_gen',aptofem_stored_keys)

  ! Set elements of mesh for cavities
    
  call set_element_cavity_no(fe_data_struct(1)%mesh_data)

! Create Solution

  ! Fluid solution
  
  call create_fe_solution(fe_data_struct(1)%soln_data,fe_data_struct(1)%mesh_data,trim(fe_solution_section), &
		aptofem_stored_keys,anal_soln)

  ! File name
    
  write(Re_str, '(es16.6)') Re
  write(DaRecipr_str, '(es16.4)') Da_reciprocal
  write(f_name_to_write, '(i0, a)') &
    id, '_Re'//trim(adjustl(Re_str))// &
    '_DaRecipr'//trim(adjustl(DaRecipr_str))//'_energy_markers.csv'

! Save markers info into CSV format
  close(unit=101)
  open(unit=101, file=f_name_to_write, status='replace')

! INFO: tX corresponds to column X, Nx corresponds to N spaces
	write(101,'("t,&
        domain_avg_v_mag,domain_avg_p,&
        total_energy_in,total_energy_out,&
        KE_in,KE_out,&
        energy_flux_in,energy_flux_out,&
        KE_flux_in,KE_flux_out,&
        p_drop,&
        face_max_p,&
        face_avg_total_energy_in,face_avg_total_energy_out,&
        face_avg_KE_in,face_avg_KE_out,&
        face_avg_energy_flux_in,face_avg_energy_flux_out,&
        face_avg_KE_flux_in,face_avg_KE_flux_out,&
        face_avg_p_drop,&
        face_avg_face_max_p,&
        all_face_avg_total_energy_in,all_face_avg_total_energy_out,&
        all_face_avg_KE_in,all_face_avg_KE_out,&
        all_face_avg_energy_flux_in,all_face_avg_energy_flux_out,&
        all_face_avg_KE_flux_in,all_face_avg_KE_flux_out,&
        all_face_avg_p_drop,&
        total_energy_diff,KE_diff,&
        total_energy_flux_diff,KE_flux_diff,&
        normalised_total_energy_diff,normalised_KE_diff,&
        normalised_total_energy_flux_diff,normalised_KE_flux_diff,&
        total_energy_frac_out_d_in,KE_frac_out_d_in,&
        total_energy_flux_frac_out_d_in,KE_flux_frac_out_d_in,&
        vol,vol_IVS,vol_slow_flow,vol_fast_flow,IVS_avg_v_mag,&
        total_viscous_dissipiation,total_viscous_dissipiation_IVS,&
        actual_viscous_dissipiation,actual_viscous_dissipiation_IVS,&
        total_rate_of_work_done,total_rate_of_work_done_IVS,&
        darcy_stress,darcy_stress_IVS,&
        corrected_drag_term,corrected_drag_term_IVS")')

	do i = 1,no_files
		
		! Fluid solution

		call read_solution_for_restart(fe_data_struct(1)%mesh_data,fe_data_struct(1)%soln_data,file_no_arr(i), &
				'soln_1',2,'../flow/output/restart/')

		! Markers

		call energy_marker(fe_data_struct(1)%mesh_data,fe_data_struct(1)%soln_data,no_energy_markers,energy_markers(:,i), &
      no_arteries,no_veins,artery_areas,vein_areas)

  end do

  do i = 1,no_energy_markers
    call calc_time_avg(time_step,no_files,energy_markers(i,1:no_files),energy_markers(i,no_files+1))
  end do

  do i = 1,no_files
		curr_t = (i-1)*time_step
    write(101,'(es16.8,a)', advance='no') curr_t,','
    do j = 1,no_energy_markers-1
      write(101,'(es16.8,a)', advance='no') energy_markers(j,i), ','
    end do
    write(101,'(es16.8)') energy_markers(no_energy_markers,i)
  end do

  ! Separator between those values and time averages
  write(101,'(a)', advance='no') '-,'
  do j = 1, no_energy_markers-1
    write(101,'(a)', advance='no') '-,'
  end do
  write(101,'(a)') '-'
    
  ! Time averages
  write(101,'(a)', advance='no') 't_avg,'
  do j = 1, no_energy_markers-1
    write(101,'(es16.8,a)', advance='no') energy_markers(j,no_files+1), ','
  end do
  write(101,'(es16.8)') energy_markers(no_energy_markers,no_files+1)
	
	close(unit=101)
	
! Delete data structures

  call delete_coupled_data(scheme_coupled_data)
  call delete_fe_data_structure(fe_data_struct,scheme_coupled_data)
  call delete_bdry_storage()
  call delete_cvty_storage()
  call delete_otlt_storage()

  deallocate(artery_radii,vein_radii,artery_areas,vein_areas)

! Finalize package

  call AptoFEM_finalize(aptofem_stored_keys)

end program flow_markers

subroutine trap_rule(h,n,f,integral)
  use param
  real(db), intent(in) :: h
  integer, intent(in) :: n
  real(db), dimension(n), intent(in) :: f
  real(db), intent(out) :: integral
  integer :: i
  integral = 0.0_db
  do i = 1,n-1
    integral = integral + 0.5_db*h*(f(i)+f(i+1))
  end do
end subroutine trap_rule

subroutine calc_time_avg(h,n,f,avg)
  use param
  real(db), intent(in) :: h
  integer, intent(in) :: n
  real(db), dimension(n), intent(in) :: f
  real(db), intent(out) :: avg
  call trap_rule(h,n,f,avg)
  avg = avg/(real(n-1,db)*h)
end subroutine calc_time_avg

!--------------------------------------------------------------------
subroutine energy_marker(mesh_data,soln_data,no_energy_markers,energy_markers, &
    no_arteries,no_veins,artery_areas,vein_areas)
!--------------------------------------------------------------------
  use param
  use fe_mesh
  use fe_solution
  use basis_fns_storage_type
  use aptofem_fe_matrix_assembly
  use problem_options

  use geometry_routines

  implicit none

  !real(db), dimension(4), intent(out) :: errors
  type(mesh), intent(inout) :: mesh_data !*FD FE mesh
  type(solution), intent(inout) :: soln_data !*FD FE solution
  integer, intent(in) :: no_energy_markers
  real(db), dimension(no_energy_markers), intent(inout) :: energy_markers

! Local variables

  type(basis_storage) :: fe_basis_info
  character(len=aptofem_length_key_def) :: control_parameter
  integer :: no_eles,no_nodes,no_faces,problem_dim,no_pdes, &
       i,j,k,qk,iv,no_quad_points,npinc, &
       no_quad_points_volume_max,no_quad_points_face_max, &
       bdry_face,dim_soln_coeff
  real(db), dimension(:,:), allocatable :: global_points_ele
  real(db), dimension(:), allocatable :: quad_weights_ele,jacobian
  real(db), dimension(:,:), allocatable :: gradient_u
  real(db), dimension(:), allocatable :: u
  real(db) :: l2_norm,h1_semi_norm,dg_norm,full_dispenal,h2_norm
  real(db), dimension(:), allocatable :: quad_weights_face,face_jacobian, &
       dispenal
  real(db), dimension(:,:), allocatable :: global_points_face,face_normals
  integer, dimension(:,:), allocatable :: global_dof_numbers1,global_dof_numbers2, &
       global_dof_numbers
  integer, dimension(2) :: loc_face_no,neighbors
  integer, dimension(:), allocatable :: no_dofs_per_variable1, &
       no_dofs_per_variable2,no_dofs_per_variable
  
  real(db) :: tol,flux,vol,vol_IVS,vol_slow_flow,vol_fast_flow, &
    total_vein_area,total_artery_area,temp_val,eff_Dm, &
    div_u,eff_Da_reciprocal,flow_spd,visc_diss,pres_work, &
    drag_term,flow_spd_square
  
  integer, intent(in) :: no_arteries,no_veins
  real(db), dimension(no_arteries), intent(in) :: artery_areas
  real(db), dimension(no_veins), intent(in) :: vein_areas

  real(db), dimension(no_arteries) :: face_pressure_in,face_pressure_face_avg_in
  real(db), dimension(no_veins) :: face_pressure_out,face_pressure_face_avg_out

  tol = 1.0d-6

  total_artery_area = sum(artery_areas)
  total_vein_area = sum(vein_areas)
  vol = 0.0_db
  vol_IVS = 0.0_db
  vol_slow_flow = 0.0_db
  vol_fast_flow = 0.0_db

  face_pressure_in = 0.0_db
  face_pressure_face_avg_in = 0.0_db
  face_pressure_out = 0.0_db
  face_pressure_face_avg_out = 0.0_db

  ! FE info
  dim_soln_coeff = get_dim_soln_coeff(soln_data)
  no_pdes = get_no_pdes(soln_data)

  call get_mesh_info(no_eles,no_nodes,no_faces,problem_dim, &
       mesh_data)

  npinc = 4
  call compute_max_no_quad_points(no_quad_points_volume_max, &
    no_quad_points_face_max,mesh_data,soln_data,npinc)

  control_parameter = 'fo_deriv_uh_ele'
  call initialize_fe_basis_storage(fe_basis_info,control_parameter,soln_data, &
    problem_dim,no_quad_points_volume_max,no_quad_points_face_max)

  allocate(gradient_u(no_pdes,problem_dim))
  allocate(u(no_pdes))
  allocate(no_dofs_per_variable(dim_soln_coeff))
  allocate(global_points_ele(problem_dim,no_quad_points_volume_max))
  allocate(jacobian(no_quad_points_volume_max))
  allocate(quad_weights_ele(no_quad_points_volume_max))
  allocate(global_dof_numbers(dim_soln_coeff,no_ele_dofs_per_var_max))

  do k = 1,no_eles

    call element_integration_info(dim_soln_coeff,problem_dim,mesh_data, &
      soln_data,k,npinc,no_quad_points_volume_max, &
      no_quad_points,global_points_ele,jacobian,quad_weights_ele, &
      global_dof_numbers,no_dofs_per_variable,fe_basis_info)

    do qk = 1,no_quad_points

      ! Used as a check whether in IVS or not
      ! Note as flow should be independent of Dm, not using actual Dm just a check of >0 or not
      call eval_tanh_smooth( &
        problem_dim,k,global_points_ele(:,qk),1.0_db,eff_Dm)
      call eval_tanh_smooth( &
        problem_dim,k,global_points_ele(:,qk),Da_reciprocal,eff_Da_reciprocal)

      u = uh_element(fe_basis_info,no_pdes,qk)
      do iv = 1,no_pdes
        gradient_u(iv,:) = grad_uh_element(fe_basis_info,problem_dim,iv,qk,1)
      end do
      div_u = 0.0_db
      do iv = 1,problem_dim
        div_u = div_u + gradient_u(iv,iv)
      end do

      flow_spd = sqrt(dot_product(u(1:problem_dim),u(1:problem_dim)))
      flow_spd_square = dot_product(u(1:problem_dim),u(1:problem_dim))

      drag_term = 2.0_db*eff_Da_reciprocal*flow_spd_square

      pres_work = u(problem_dim+1)*div_u
      do i = 1,problem_dim
        do j = 1,problem_dim
          visc_diss = gradient_u(i,j)*gradient_u(i,j) &
            +gradient_u(i,j)*gradient_u(j,i) &
            +gradient_u(j,i)*gradient_u(i,j) &
            +gradient_u(j,i)*gradient_u(j,i)
        end do
      end do

      energy_markers(1) = energy_markers(1) + &
        jacobian(qk)*quad_weights_ele(qk)*flow_spd

      ! (1/Da)*|u|
      energy_markers(55) = energy_markers(55) + &
        jacobian(qk)*quad_weights_ele(qk)*eff_Da_reciprocal*flow_spd

      ! (2/Da)*|u|^2
      energy_markers(57) = energy_markers(57) + &
        jacobian(qk)*quad_weights_ele(qk)*drag_term

      energy_markers(2) = energy_markers(2) + &
        jacobian(qk)*quad_weights_ele(qk)* &
          u(problem_dim+1)

      energy_markers(49) = energy_markers(49) - &
        jacobian(qk)*quad_weights_ele(qk)*pres_work
      do i = 1,problem_dim
        do j = 1,problem_dim
          energy_markers(49) = energy_markers(49) + &
            jacobian(qk)*quad_weights_ele(qk)* &
              (gradient_u(i,j)*gradient_u(j,i)+gradient_u(j,i)*gradient_u(j,i))
        end do
      end do
      energy_markers(51) = energy_markers(51) + &
        jacobian(qk)*quad_weights_ele(qk)* &
          0.5_db*visc_diss
      energy_markers(53) = energy_markers(53) + &
        jacobian(qk)*quad_weights_ele(qk)* &
          (0.5_db*visc_diss - pres_work)

      vol = vol + jacobian(qk)*quad_weights_ele(qk)

      if (eff_Dm > tol) then
        vol_IVS = vol_IVS + jacobian(qk)*quad_weights_ele(qk)

        energy_markers(48) = energy_markers(48) + &
          jacobian(qk)*quad_weights_ele(qk)*flow_spd

        ! (1/Da)*|u|
        energy_markers(56) = energy_markers(56) + &
          jacobian(qk)*quad_weights_ele(qk)*eff_Da_reciprocal*flow_spd

        ! (2/Da)*|u|^2
        energy_markers(58) = energy_markers(58) + &
          jacobian(qk)*quad_weights_ele(qk)*drag_term

        energy_markers(50) = energy_markers(50) - &
          jacobian(qk)*quad_weights_ele(qk)*u(problem_dim+1)*div_u
        do i = 1,problem_dim
          do j = 1,problem_dim
            energy_markers(50) = energy_markers(50) + &
              jacobian(qk)*quad_weights_ele(qk)* &
                (gradient_u(i,j)*gradient_u(j,i)+gradient_u(j,i)*gradient_u(j,i))
          end do
        end do
        energy_markers(52) = energy_markers(52) + &
          jacobian(qk)*quad_weights_ele(qk)* &
            0.5_db*visc_diss
        energy_markers(54) = energy_markers(54) + &
          jacobian(qk)*quad_weights_ele(qk)* &
            (0.5_db*visc_diss - pres_work)

      end if

        ! This vel from Dellschaft 2020: slow flow defined 5x10^-4 m/s, V_c = 1.5x10^-1 m/s, fast 1x10^-3 m/s
      if (flow_spd < (1.0_db / 300.0_db)) then
        vol_slow_flow = vol_slow_flow + jacobian(qk)*quad_weights_ele(qk)
      else if (flow_spd > (1.0_db / 150.0_db)) then
        vol_fast_flow = vol_fast_flow + jacobian(qk)*quad_weights_ele(qk)
      end if

     end do

  end do

  energy_markers(44) = vol
  energy_markers(45) = vol_IVS
  energy_markers(46) = vol_slow_flow
  energy_markers(47) = vol_fast_flow

  ! Avg. IVS flow
  energy_markers(48) = energy_markers(48) / vol_IVS
  ! Avg. total flow, pressure
  energy_markers(1:2) = energy_markers(1:2) / vol

  call delete_fe_basis_storage(fe_basis_info)

  deallocate(gradient_u,no_dofs_per_variable,global_points_ele,jacobian, &
    quad_weights_ele,global_dof_numbers)



! -------------------- FACE --------------------

  control_parameter = 'uh_face'
  call initialize_fe_basis_storage(fe_basis_info,control_parameter,soln_data, &
       problem_dim,no_quad_points_volume_max,no_quad_points_face_max)

  allocate(global_points_face(problem_dim,no_quad_points_face_max))
  allocate(face_jacobian(no_quad_points_face_max))
  allocate(face_normals(problem_dim,no_quad_points_face_max))
  allocate(quad_weights_face(no_quad_points_face_max))
  allocate(no_dofs_per_variable1(dim_soln_coeff))
  allocate(no_dofs_per_variable2(dim_soln_coeff))
  allocate(global_dof_numbers1(dim_soln_coeff,no_ele_dofs_per_var_max))
  allocate(global_dof_numbers2(dim_soln_coeff,no_ele_dofs_per_var_max))

  do k = 1,no_faces
        
    call face_integration_info(dim_soln_coeff,problem_dim,mesh_data,soln_data, &
         k,neighbors,loc_face_no,npinc,no_quad_points_face_max, &
         no_quad_points,global_points_face,face_jacobian,face_normals, &
         quad_weights_face,global_dof_numbers1,no_dofs_per_variable1, &
         bdry_face,global_dof_numbers2,no_dofs_per_variable2, &
         fe_basis_info)

! Boundary face
	 
! 1: (1/A) int_{out} (p + 0.5*v^2)*dot(v,n), (1/A) int_{in} (p + 0.5*v^2)*dot(v,n) [energy_flux]

	 if (bdry_face > 1 .AND. bdry_face <= 100) then

        do qk = 1,no_quad_points
		
		      u = uh_face1(fe_basis_info,no_pdes,qk)
		   
		      flux = abs(dot_product(u(1:problem_dim),face_normals(1:problem_dim,qk)))
          
          temp_val = face_jacobian(qk)*quad_weights_face(qk)* &
            ((1.0_db/Re)*u(no_pdes) + &
              0.5_db*dot_product(u(1:problem_dim),u(1:problem_dim)))
          
          energy_markers(3) = energy_markers(3) + &
            temp_val

          energy_markers(7) = energy_markers(7) + &
            temp_val*flux

          energy_markers(13) = energy_markers(13) + &
            (1.0_db/artery_areas(bdry_face-1))*temp_val

          energy_markers(17) = energy_markers(17) + &
            (1.0_db/artery_areas(bdry_face-1))*temp_val*flux

          energy_markers(23) = energy_markers(23) + &
            (1.0_db/total_artery_area)*temp_val

          energy_markers(27) = energy_markers(27) + &
            (1.0_db/total_artery_area)*temp_val*flux

        end do
		
	 end if

	 if (bdry_face > 100 .AND. bdry_face <= 200) then

    do qk = 1,no_quad_points

      u = uh_face1(fe_basis_info,no_pdes,qk)
   
      flux = abs(dot_product(u(1:problem_dim),face_normals(1:problem_dim,qk)))
          
      temp_val = face_jacobian(qk)*quad_weights_face(qk)* &
        ((1.0_db/Re)*u(no_pdes) + &
          0.5_db*dot_product(u(1:problem_dim),u(1:problem_dim)))
      
      energy_markers(4) = energy_markers(4) + &
        temp_val

      energy_markers(8) = energy_markers(8) + &
        temp_val*flux

      energy_markers(14) = energy_markers(14) + &
        (1.0_db/vein_areas(bdry_face-100))*temp_val

      energy_markers(18) = energy_markers(18) + &
        (1.0_db/vein_areas(bdry_face-100))*temp_val*flux

      energy_markers(24) = energy_markers(24) + &
        (1.0_db/total_vein_area)*temp_val

      energy_markers(28) = energy_markers(28) + &
        (1.0_db/total_vein_area)*temp_val*flux

    end do

   end if
	 
! 2: 0.5*v^2*dot(v,n) [KE flux]

	 if (bdry_face > 1 .AND. bdry_face <= 100) then

      do qk = 1,no_quad_points

        u = uh_face1(fe_basis_info,no_pdes,qk)
    
        flux = abs(dot_product(u(1:problem_dim),face_normals(1:problem_dim,qk)))
        
        temp_val = face_jacobian(qk)*quad_weights_face(qk)* &
          0.5_db*dot_product(u(1:problem_dim),u(1:problem_dim))
        
        energy_markers(5) = energy_markers(5) + &
          temp_val

        energy_markers(9) = energy_markers(9) + &
          temp_val*flux

        energy_markers(15) = energy_markers(15) + &
          (1.0_db/artery_areas(bdry_face-1))*temp_val

        energy_markers(19) = energy_markers(19) + &
          (1.0_db/artery_areas(bdry_face-1))*temp_val*flux

        energy_markers(25) = energy_markers(25) + &
          (1.0_db/total_artery_area)*temp_val

        energy_markers(29) = energy_markers(29) + &
          (1.0_db/total_artery_area)*temp_val*flux

      end do

    end if

    if (bdry_face > 100 .AND. bdry_face <= 200) then

      do qk = 1,no_quad_points

        u = uh_face1(fe_basis_info,no_pdes,qk)

        flux = abs(dot_product(u(1:problem_dim),face_normals(1:problem_dim,qk)))
            
        temp_val = face_jacobian(qk)*quad_weights_face(qk)* &
          0.5_db*dot_product(u(1:problem_dim),u(1:problem_dim))
        
        energy_markers(6) = energy_markers(6) + &
          temp_val

        energy_markers(10) = energy_markers(10) + &
          temp_val*flux

        energy_markers(16) = energy_markers(16) + &
          (1.0_db/vein_areas(bdry_face-100))*temp_val

        energy_markers(20) = energy_markers(20) + &
          (1.0_db/vein_areas(bdry_face-100))*temp_val*flux

        energy_markers(26) = energy_markers(26) + &
          (1.0_db/total_vein_area)*temp_val

        energy_markers(30) = energy_markers(30) + &
          (1.0_db/total_vein_area)*temp_val*flux

      end do

    end if

    ! 3: Pressure

    if (bdry_face > 1 .AND. bdry_face <= 100) then

      do qk = 1,no_quad_points

        u = uh_face1(fe_basis_info,no_pdes,qk)
        
        temp_val = face_jacobian(qk)*quad_weights_face(qk)*u(no_pdes)
        
        energy_markers(11) = energy_markers(11) + &
          temp_val

        ! For marker 12, face_max_p
        face_pressure_in(bdry_face-1) = face_pressure_in(bdry_face-1) + &
          temp_val

        ! face_avg_p_drop
        energy_markers(21) = energy_markers(21) + &
          (1.0_db/artery_areas(bdry_face-1))*temp_val

        ! For marker 22, face_avg_face_max_p
        face_pressure_face_avg_in(bdry_face-1) = face_pressure_face_avg_in(bdry_face-1) + &
          (1.0_db/artery_areas(bdry_face-1))*temp_val

        ! all_face_avg_p_drop
        energy_markers(31) = energy_markers(31) + &
          (1.0_db/total_artery_area)*temp_val

      end do

    end if

    if (bdry_face > 100 .AND. bdry_face <= 200) then

      do qk = 1,no_quad_points

        u = uh_face1(fe_basis_info,no_pdes,qk)
          
        temp_val = face_jacobian(qk)*quad_weights_face(qk)*u(no_pdes)
        
        energy_markers(11) = energy_markers(11) - &
          temp_val

        ! For marker 12, face_max_p
        face_pressure_out(bdry_face-100) = face_pressure_out(bdry_face-100) + &
          temp_val

        ! face_avg_p_drop
        energy_markers(21) = energy_markers(21) - &
          (1.0_db/vein_areas(bdry_face-100))*temp_val

        ! For marker 22, face_avg_face_max_p
        face_pressure_face_avg_out(bdry_face-100) = face_pressure_face_avg_out(bdry_face-100) + &
          (1.0_db/vein_areas(bdry_face-100))*temp_val

        ! all_face_avg_p_drop
        energy_markers(31) = energy_markers(31) - &
          (1.0_db/total_vein_area)*temp_val

      end do

    end if

  end do
  
  deallocate(global_points_face,face_jacobian,face_normals,quad_weights_face, &
       no_dofs_per_variable1,no_dofs_per_variable2, &
       global_dof_numbers1,global_dof_numbers2,u)

  call delete_fe_basis_storage(fe_basis_info)
  
  ! face_max_p, greatest p diff
  energy_markers(12) = maxval(face_pressure_in)-minval(face_pressure_out)
  ! face_avg_face_max_p, average of the greatest p diff
  energy_markers(22) = maxval(face_pressure_face_avg_in)-minval(face_pressure_face_avg_out)

  ! total_energy_diff
  energy_markers(32) = energy_markers(3)-energy_markers(4)
  ! KE_diff
  energy_markers(33) = energy_markers(5)-energy_markers(6)
  ! total_energy_flux_diff
  energy_markers(34) = energy_markers(7)-energy_markers(8)
  ! KE_flux_diff
  energy_markers(35) = energy_markers(9)-energy_markers(10)
  ! normalised_total_energy_diff
  energy_markers(36) = energy_markers(32)/energy_markers(3)
  ! normalised_KE_diff
  energy_markers(37) = energy_markers(33)/energy_markers(5)
  ! normalised_total_energy_flux_diff
  energy_markers(38) = energy_markers(34)/energy_markers(7)
  ! normalised_KE_flux_diff
  energy_markers(39) = energy_markers(35)/energy_markers(9)
  ! total_energy_frac_out_d_in
  energy_markers(40) = energy_markers(4)/energy_markers(3)
  ! KE_frac_out_d_in
  energy_markers(41) = energy_markers(6)/energy_markers(5)
  ! total_energy_flux_frac_out_d_in
  energy_markers(42) = energy_markers(8)/energy_markers(7)
  ! KE_flux_frac_out_d_in
  energy_markers(43) = energy_markers(10)/energy_markers(9)

  do i = 1,no_energy_markers
    write(*,fmt='(A,i0,A,ES15.3E3)') "energy_markers(",i,") = ",energy_markers(i)
  end do

end subroutine energy_marker
