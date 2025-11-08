! To dimensionalise: replace Re with rho, 1.000_db with mu
module navier_stokes_jac_matrix_residual

  use coupled_data_type

  use general_routines
  use geometry_routines

contains

!--------------------------------------------------------------------
  subroutine ns_element_nonlin_residual(element_rhs, &
                                        mesh_data, soln_data, facet_data, fe_basis_info)
!--------------------------------------------------------------------
    use problem_options
    use coupled_data_storage

    include 'assemble_residual_element.h'

! Local variables

    integer :: qk, i, j, ieqn
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points) &
      :: floc
    real(db), dimension(facet_data%problem_dim, facet_data%problem_dim, &
                        facet_data%no_quad_points) :: fluxes
!  real(db) :: div_u
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points) :: &
      interpolant_uh
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points, &
                        facet_data%problem_dim) :: gradient_uh
    real(db), dimension(facet_data%dim_soln_coeff, facet_data%no_quad_points, &
                        facet_data%problem_dim, maxval(facet_data%no_dofs_per_variable)) :: grad_phi
    real(db), dimension(facet_data%dim_soln_coeff, facet_data%no_quad_points, &
                        maxval(facet_data%no_dofs_per_variable)) :: phi
    real(db), dimension(facet_data%no_quad_points) :: div_u, eff_permeability_reciprocal

    integer :: prev_sol, bdf_order
    real(db) :: bdf_scaling_factor, current_time, time_step
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points, &
                        facet_data%scheme_user_data%bdf_order) :: uh_previous_time_step
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points) :: &
      soln_extrap
    real(db), dimension(facet_data%problem_dim, facet_data%no_quad_points) :: &
      vel_bdf_source

    associate ( &
      dim_soln_coeff => facet_data%dim_soln_coeff, &
      no_pdes => facet_data%no_pdes, &
      problem_dim => facet_data%problem_dim, &
      no_quad_points => facet_data%no_quad_points, &
      global_points_ele => facet_data%global_points, &
      integral_weighting => facet_data%integral_weighting, &
      element_number => facet_data%element_number, &
      element_region_id => facet_data%element_region_id, &
      no_dofs_per_variable => facet_data%no_dofs_per_variable, &
      global_dof_numbers => facet_data%global_dof_numbers, &
      scheme_user_data => facet_data%scheme_user_data)

      bdf_order = scheme_user_data%bdf_order
      bdf_scaling_factor = scheme_user_data%bdf_scaling_factor
      current_time = scheme_user_data%current_time
      time_step = scheme_user_data%time_step

      select type (scheme_user_data)
      type is (coupled_user_data)
        do prev_sol = 1, bdf_order
          call compute_uh_with_basis_fns_pts_from_array(uh_previous_time_step(:, :, prev_sol), &
                                                        no_pdes, no_quad_points, dim_soln_coeff, no_dofs_per_variable, &
                                                        global_dof_numbers, fe_basis_info%basis_element, soln_data, &
                                                        soln_data%no_dofs, fe_data_struct(1)%soln_data_prev(prev_sol)%soln_values)
        end do
      class default
        do prev_sol = 1, bdf_order
          call compute_uh_with_basis_fns_pts_from_array(uh_previous_time_step(:, :, prev_sol), &
                                                        no_pdes, no_quad_points, dim_soln_coeff, no_dofs_per_variable, &
                                                        global_dof_numbers, fe_basis_info%basis_element, soln_data, &
                                                        soln_data%no_dofs, scheme_user_data%uh_previous_time_step_ms(:, prev_sol))
        end do
      end select

      element_rhs = 0.0_db

      do qk = 1, no_quad_points

        interpolant_uh(1:no_pdes, qk) = uh_element(fe_basis_info, no_pdes, qk)
        do i = 1, no_pdes
          gradient_uh(i, qk, 1:problem_dim) = grad_uh_element(fe_basis_info, problem_dim, i, qk, 1)
        end do
        call forcing_fn(floc(:, qk), global_points_ele(:, qk), problem_dim, no_pdes, current_time)
        call bdf_extrapolate(soln_extrap(:, qk), uh_previous_time_step(:, qk, :), no_pdes, bdf_order)
        call bdf_source(vel_bdf_source(:, qk), uh_previous_time_step(1:problem_dim, qk, :), &
                        problem_dim, bdf_order, time_step)
        call extrap_adv_fluxes(interpolant_uh(:, qk), fluxes(:, :, qk), problem_dim, no_pdes, &
                               soln_extrap(:, qk))
        call eval_tanh_smooth(problem_dim, element_number, global_points_ele(:, qk), &
                              Da_reciprocal, eff_permeability_reciprocal(qk))
        div_u(qk) = 0.0_db

        do ieqn = 1, problem_dim
          div_u(qk) = div_u(qk) + gradient_uh(ieqn, qk, ieqn)
        end do

      end do

! Calculate Basis Functions

      do i = 1, dim_soln_coeff
        grad_phi(i, 1:no_quad_points, 1:problem_dim, 1:no_dofs_per_variable(i)) = fe_basis_info%basis_element &
                                        %deriv_basis_fns(i)%grad_data(1:no_quad_points, 1:problem_dim, 1:no_dofs_per_variable(i), 1)
        phi(i, 1:no_quad_points, 1:no_dofs_per_variable(i)) = fe_basis_info%basis_element%basis_fns(i) &
                                                              %fem_basis_fns(1:no_quad_points, 1:no_dofs_per_variable(i), 1)
      end do

! Momentum Equations

      do ieqn = 1, problem_dim

! Loop over quadrature points

        do qk = 1, no_quad_points

! Loop over phi_i

          do i = 1, no_dofs_per_variable(ieqn)

            ! First: time terms, (-) for RHS; Second: grad v terms (adv, pressure);
            ! Third: v terms (forcing !!!needs Froude no.!!!); Fourth: BDF vel source
            ! Fifth: udiv(u) stab
            element_rhs(ieqn, i) = element_rhs(ieqn, i) &
                                   + integral_weighting(qk)*( &
                                   Re*( &
                                   (-bdf_scaling_factor*interpolant_uh(ieqn, qk) &
                                    + vel_bdf_source(ieqn, qk))*phi(ieqn, qk, i) &
                                   + dot_product(fluxes(1:problem_dim, ieqn, qk), &
                                                 grad_phi(ieqn, qk, 1:problem_dim, i)) &
                                   ) &
                        - 0.5_db*1.000_db*cal_sym_stress(gradient_uh(:, qk, :), grad_phi(:, qk, :, i), ieqn, problem_dim, no_pdes) &
                                   + interpolant_uh(problem_dim + 1, qk)*grad_phi(ieqn, qk, ieqn, i) &
                                   + floc(ieqn, qk)*phi(ieqn, qk, i))

            !+ div_u(qk)*interpolant_uh(ieqn,qk)*phi(ieqn,qk,i)/2.0_db)

            ! Darcy term
            element_rhs(ieqn, i) = element_rhs(ieqn, i) &
                                   + integral_weighting(qk)*( &
                                   -1.000_db*eff_permeability_reciprocal(qk)*interpolant_uh(ieqn, qk)*phi(ieqn, qk, i) &
                                   )

          end do
        end do
      end do

! Divergence-free constraint

! Loop over quadrature points

      do qk = 1, no_quad_points

! Loop over phi_i

        do i = 1, no_dofs_per_variable(problem_dim + 1)

          element_rhs(no_pdes, i) = element_rhs(no_pdes, i) + integral_weighting(qk) &
                                    *(floc(problem_dim + 1, qk) - div_u(qk))*phi(problem_dim + 1, qk, i)
        end do
      end do

    end associate

  end subroutine ns_element_nonlin_residual

!--------------------------------------------------------------------
  subroutine ns_dg_face_nonlin_residual(face_residual_p, face_residual_m, &
                                        mesh_data, soln_data, facet_data, fe_basis_info)
!--------------------------------------------------------------------
    use problem_options
    use coupled_data_storage
    use aptofem_fe_matrix_assembly

    include 'assemble_residual_int_bdry_face.h'

! Local variables

    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points) &
      :: uloc
    real(db), dimension(facet_data%no_quad_points) &
      :: p_bc
    real(db), dimension(facet_data%problem_dim, &
                        facet_data%no_quad_points) :: unloc, nflxsoln
    integer :: i, qk, ieqn
    real(db) :: full_dispenal
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points) :: &
      interpolant_uh1, interpolant_uh2
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points, &
                        facet_data%problem_dim) :: gradient_uh1, gradient_uh2
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points, &
                        facet_data%problem_dim, maxval(facet_data%no_dofs_per_variable1)) &
      :: grad_phi1
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points, &
                        facet_data%problem_dim, maxval(facet_data%no_dofs_per_variable2)) &
      :: grad_phi2
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points, &
                        maxval(facet_data%no_dofs_per_variable1)) :: phi1
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points, &
                        maxval(facet_data%no_dofs_per_variable2)) :: phi2
    real(db), dimension(facet_data%dim_soln_coeff) :: dispenal_local
    real(db) :: lf_flux
    logical :: bdf_semi_implicit

    integer :: prev_sol, bdf_order
    real(db) :: current_time, transmission_stress, time_step, udotn_factor_bdf
    real(db), dimension(facet_data%no_quad_points) :: bf_stab_term, total_normal_stress
    real(db), dimension(facet_data%no_quad_points, &
                        facet_data%scheme_user_data%bdf_order) :: normal_stress
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points, &
                        facet_data%scheme_user_data%bdf_order) :: uh_previous_time_step1
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points, &
                        facet_data%scheme_user_data%bdf_order) :: uh_previous_time_step2
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points, &
                        facet_data%problem_dim, facet_data%scheme_user_data%bdf_order) :: &
      gradient_uh_previous_time_step1
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points, &
                        facet_data%problem_dim, facet_data%scheme_user_data%bdf_order) :: &
      gradient_uh_previous_time_step2
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points) :: &
      soln_extrap1
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points) :: &
      soln_extrap2
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points) :: &
      uh_previous_iterate1
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points, &
                        facet_data%problem_dim) :: gradient_uh_previous_iterate1
    real(db), dimension(facet_data%problem_dim, facet_data%no_quad_points) :: &
      interpolant_vel1_normal, interpolant_vel_extrap1_normal
    real(db), dimension(facet_data%problem_dim, facet_data%no_quad_points) :: &
      interpolant_vel1_normal_prev !!!!ONLY BDF1 - TEMPORARY
    real(db), dimension(facet_data%problem_dim, facet_data%no_quad_points) :: &
      tang_stress_vector_u
    real(db), dimension(facet_data%problem_dim, facet_data%problem_dim, &
                        facet_data%no_quad_points) :: stress_tensor, surface_tensor
    real(db), dimension(facet_data%problem_dim, facet_data%no_quad_points, &
                        maxval(facet_data%no_dofs_per_variable1)) :: tang_penalty, tang_stress_vector_phi

    associate ( &
      dim_soln_coeff => facet_data%dim_soln_coeff, &
      no_pdes => facet_data%no_pdes, &
      problem_dim => facet_data%problem_dim, &
      no_quad_points => facet_data%no_quad_points, &
      global_points_face => facet_data%global_points, &
      integral_weighting => facet_data%integral_weighting, &
      face_number => facet_data%face_number, &
      neighbours => facet_data%neighbours, &
      interior_face_boundary_no => facet_data%interior_face_boundary_no, &
      face_element_region_ids => facet_data%face_element_region_ids, &
      bdry_face => facet_data%bdry_no, &
      no_dofs_per_variable1 => facet_data%no_dofs_per_variable1, &
      no_dofs_per_variable2 => facet_data%no_dofs_per_variable2, &
      face_normals => facet_data%face_normals, &
      dispenal => facet_data%dispenal, &
      global_dof_numbers1 => facet_data%global_dof_numbers1, &
      global_dof_numbers2 => facet_data%global_dof_numbers2, &
      scheme_user_data => facet_data%scheme_user_data)

      bdf_order = scheme_user_data%bdf_order
      current_time = scheme_user_data%current_time
      bdf_semi_implicit = scheme_user_data%bdf_semi_implicit
      time_step = scheme_user_data%time_step

      face_residual_p = 0.0_db
      face_residual_m = 0.0_db

      dispenal_local = dispenal

      full_dispenal = 1.000_db*interior_penalty_parameter*dispenal_local(1)

      if (bdry_face > 0) then

        select type (scheme_user_data)
        type is (coupled_user_data)
          do prev_sol = 1, bdf_order
            if (bdry_face > 300 .AND. bdry_face <= 400) then
              call compute_uh_gradient_uh_with_basis_fns_pts_from_array( &
                uh_previous_time_step1(:, :, prev_sol), gradient_uh_previous_time_step1(:, :, :, prev_sol), &
                problem_dim, no_pdes, no_quad_points, dim_soln_coeff, facet_data%no_dofs_per_variable1, &
                facet_data%global_dof_numbers1, fe_basis_info%basis_face1, soln_data, &
                soln_data%no_dofs, fe_data_struct(1)%soln_data_prev(prev_sol)%soln_values)
            else
              call compute_uh_with_basis_fns_pts_from_array(uh_previous_time_step1(:, :, prev_sol), &
                                                        no_pdes, no_quad_points, dim_soln_coeff, facet_data%no_dofs_per_variable1, &
                                                            facet_data%global_dof_numbers1, fe_basis_info%basis_face1, soln_data, &
                                                          soln_data%no_dofs, fe_data_struct(1)%soln_data_prev(prev_sol)%soln_values)
            end if
          end do
        class default
          if (bdry_face > 300) then
            write (io_err, *) 'ERROR: bdry_face > 300 can only be used with coupled_user_data type'
            stop
          end if
          do prev_sol = 1, bdf_order
            call compute_uh_with_basis_fns_pts_from_array(uh_previous_time_step1(:, :, prev_sol), &
                                                        no_pdes, no_quad_points, dim_soln_coeff, facet_data%no_dofs_per_variable1, &
                                                          facet_data%global_dof_numbers1, fe_basis_info%basis_face1, soln_data, &
                                                          soln_data%no_dofs, scheme_user_data%uh_previous_time_step_ms(:, prev_sol))
          end do
        end select

! Calculate value of solution at quadrature points

        do qk = 1, no_quad_points
          interpolant_uh1(1:no_pdes, qk) = uh_face1(fe_basis_info, no_pdes, qk)
          do i = 1, no_pdes
            gradient_uh1(i, qk, 1:problem_dim) = grad_uh_face1(fe_basis_info, problem_dim, i, qk, 1)
          end do

          call compute_boundary_condition_bdf(interpolant_uh2(:, qk), &
                                              uloc(:, qk), interpolant_uh1(:, qk), face_normals(:, qk), &
                                              problem_dim, no_pdes, abs(bdry_face), global_points_face(:, qk), &
                                              scheme_user_data%current_time, unloc(:, qk), p_bc(qk))
          call lax_friedrichs_bdf_bdry(nflxsoln(:, qk), interpolant_uh1(:, qk), &
                                       interpolant_uh2(:, qk), face_normals(:, qk), problem_dim, no_pdes, bdry_face, &
                                       bdf_order, bdf_semi_implicit, uh_previous_time_step1(:, qk, :))

          call bdf_extrapolate(soln_extrap1(:, qk), &
                               uh_previous_time_step1(:, qk, :), no_pdes, bdf_order)
          call normal_vel_stab(bf_stab_term(qk), &
                               soln_extrap1(:, qk), face_normals(:, qk), &
                               problem_dim, no_pdes, backflow_stab_parameter)
        end do

        do i = 1, no_pdes
          grad_phi1(i, 1:no_quad_points, 1:problem_dim, 1:no_dofs_per_variable1(i)) = fe_basis_info%basis_face1 &
                                                   %deriv_basis_fns(i)%grad_data(1:no_quad_points, :, 1:no_dofs_per_variable1(i), 1)
          phi1(i, 1:no_quad_points, 1:no_dofs_per_variable1(i)) = fe_basis_info%basis_face1%basis_fns(i) &
                                                                  %fem_basis_fns(1:no_quad_points, 1:no_dofs_per_variable1(i), 1)
        end do

        if (abs(bdry_face) >= 1 .and. abs(bdry_face) <= 100) then

! Dirichlet boundary

! Loop over quadrature points

! =======================================================================================
! Momentum Equations

          do ieqn = 1, problem_dim
            do qk = 1, no_quad_points

! Loop over phi_i

              do i = 1, no_dofs_per_variable1(ieqn)

                face_residual_p(ieqn, i) = face_residual_p(ieqn, i) &
                                           + integral_weighting(qk)*( &
                                           Re*(-nflxsoln(ieqn, qk)*phi1(ieqn, qk, i)) &
                                           + (1.000_db*cal_bdry_sym_stress_u(gradient_uh1(:, qk, :), &
                                                                             face_normals(:, qk), ieqn, problem_dim, no_pdes) &
                                              - full_dispenal*(interpolant_uh1(ieqn, qk) - uloc(ieqn, qk))) &
                                           *phi1(ieqn, qk, i) &
                                           + 1.000_db*cal_bdry_sym_stress_v(grad_phi1(ieqn, qk, :, i), interpolant_uh1(:, qk), &
                                                                     uloc(:, qk), face_normals(:, qk), ieqn, problem_dim, no_pdes) &
                                           - interpolant_uh1(problem_dim + 1, qk)*face_normals(ieqn, qk) &
                                           *phi1(ieqn, qk, i) &
                                           )

              end do
            end do
          end do

! =======================================================================================
! Divergence free constraint

          do qk = 1, no_quad_points

! Loop over phi_i

            do i = 1, no_dofs_per_variable1(problem_dim + 1)

              face_residual_p(no_pdes, i) = face_residual_p(no_pdes, i) + integral_weighting(qk)* &
                                            dot_product(interpolant_uh1(1:problem_dim, qk) - uloc(1:problem_dim, qk), &
                                                        face_normals(:, qk))*phi1(problem_dim + 1, qk, i)

            end do
          end do

! =======================================================================================

        else if (abs(bdry_face) >= 101 .and. abs(bdry_face) <= 200) then

! Neumann Boundary

! Loop over quadrature points

! =======================================================================================
! Momentum Equations

          do ieqn = 1, problem_dim
            do qk = 1, no_quad_points

              lf_flux = interpolant_uh1(ieqn, qk)*dot_product(soln_extrap1(1:problem_dim, qk), face_normals(:, qk))

! Loop over phi_i

              do i = 1, no_dofs_per_variable1(ieqn)

                face_residual_p(ieqn, i) = face_residual_p(ieqn, i) &
                                           - integral_weighting(qk)*( &
                                           Re*(nflxsoln(ieqn, qk) + bf_stab_term(qk)*interpolant_uh1(ieqn, qk)) &
                                           + unloc(ieqn, qk) &
                                           + p_bc(qk)*face_normals(ieqn, qk) &
                                           )*phi1(ieqn, qk, i)

              end do
            end do
          end do

        end if

      else

! Interior face

        select type (scheme_user_data)
        type is (coupled_user_data)
          do prev_sol = 1, bdf_order
            call compute_uh_with_basis_fns_pts_from_array(uh_previous_time_step1(:, :, prev_sol), &
                                                        no_pdes, no_quad_points, dim_soln_coeff, facet_data%no_dofs_per_variable1, &
                                                          facet_data%global_dof_numbers1, fe_basis_info%basis_face1, soln_data, &
                                                          soln_data%no_dofs, fe_data_struct(1)%soln_data_prev(prev_sol)%soln_values)
            call compute_uh_with_basis_fns_pts_from_array(uh_previous_time_step2(:, :, prev_sol), &
                                                        no_pdes, no_quad_points, dim_soln_coeff, facet_data%no_dofs_per_variable2, &
                                                          facet_data%global_dof_numbers2, fe_basis_info%basis_face2, soln_data, &
                                                          soln_data%no_dofs, fe_data_struct(1)%soln_data_prev(prev_sol)%soln_values)
          end do
        class default
          do prev_sol = 1, bdf_order
            call compute_uh_with_basis_fns_pts_from_array(uh_previous_time_step1(:, :, prev_sol), &
                                                        no_pdes, no_quad_points, dim_soln_coeff, facet_data%no_dofs_per_variable1, &
                                                          facet_data%global_dof_numbers1, fe_basis_info%basis_face1, soln_data, &
                                                          soln_data%no_dofs, scheme_user_data%uh_previous_time_step_ms(:, prev_sol))
            call compute_uh_with_basis_fns_pts_from_array(uh_previous_time_step2(:, :, prev_sol), &
                                                        no_pdes, no_quad_points, dim_soln_coeff, facet_data%no_dofs_per_variable2, &
                                                          facet_data%global_dof_numbers2, fe_basis_info%basis_face2, soln_data, &
                                                          soln_data%no_dofs, scheme_user_data%uh_previous_time_step_ms(:, prev_sol))
          end do
        end select

        do qk = 1, no_quad_points
          interpolant_uh1(:, qk) = uh_face1(fe_basis_info, no_pdes, qk)
          interpolant_uh2(:, qk) = uh_face2(fe_basis_info, no_pdes, qk)

          do i = 1, no_pdes
            gradient_uh1(i, qk, 1:problem_dim) = grad_uh_face1(fe_basis_info, problem_dim, i, qk, 1)
            gradient_uh2(i, qk, 1:problem_dim) = grad_uh_face2(fe_basis_info, problem_dim, i, qk, 1)
          end do

          call lax_friedrichs_bdf_int(nflxsoln(:, qk), interpolant_uh1(:, qk), &
                                      interpolant_uh2(:, qk), face_normals(:, qk), problem_dim, no_pdes, bdry_face, &
                                      bdf_order, bdf_semi_implicit, uh_previous_time_step1(:, qk, :), &
                                      uh_previous_time_step2(:, qk, :))

        end do

        do i = 1, no_pdes
          grad_phi1(i, 1:no_quad_points, 1:problem_dim, 1:no_dofs_per_variable1(i)) = fe_basis_info%basis_face1 &
                                                   %deriv_basis_fns(i)%grad_data(1:no_quad_points, :, 1:no_dofs_per_variable1(i), 1)
          grad_phi2(i, 1:no_quad_points, 1:problem_dim, 1:no_dofs_per_variable2(i)) = fe_basis_info%basis_face2 &
                                                   %deriv_basis_fns(i)%grad_data(1:no_quad_points, :, 1:no_dofs_per_variable2(i), 1)
          phi1(i, 1:no_quad_points, 1:no_dofs_per_variable1(i)) = fe_basis_info%basis_face1%basis_fns(i) &
                                                                  %fem_basis_fns(1:no_quad_points, 1:no_dofs_per_variable1(i), 1)
          phi2(i, 1:no_quad_points, 1:no_dofs_per_variable2(i)) = fe_basis_info%basis_face2%basis_fns(i) &
                                                                  %fem_basis_fns(1:no_quad_points, 1:no_dofs_per_variable2(i), 1)
        end do

! Loop over quadrature points

! =======================================================================================
! Momentum Equations

        do ieqn = 1, problem_dim

          do qk = 1, no_quad_points

! Loop over phi_i

            do i = 1, no_dofs_per_variable1(ieqn)

              face_residual_p(ieqn, i) = face_residual_p(ieqn, i) + integral_weighting(qk)*( &
                                         Re*(-nflxsoln(ieqn, qk)*phi1(ieqn, qk, i)) &
                + (0.5_db*1.000_db*(cal_bdry_sym_stress_u(gradient_uh1(:, qk, :), face_normals(:, qk), ieqn, problem_dim, no_pdes) &
                                  + cal_bdry_sym_stress_u(gradient_uh2(:, qk, :), face_normals(:, qk), ieqn, problem_dim, no_pdes) &
                                                             ) &
                                            - full_dispenal*(interpolant_uh1(ieqn, qk) - interpolant_uh2(ieqn, qk))) &
                                         *phi1(ieqn, qk, i) &
                                        + 0.5_db*1.000_db*cal_bdry_sym_stress_v(grad_phi1(ieqn, qk, :, i), interpolant_uh1(:, qk), &
                                                          interpolant_uh2(:, qk), face_normals(:, qk), ieqn, problem_dim, no_pdes) &
                                         - 0.5_db*(interpolant_uh1(problem_dim + 1, qk) + interpolant_uh2(problem_dim + 1, qk))* &
                                         face_normals(ieqn, qk)*phi1(ieqn, qk, i) &
                                         )

            end do

! Loop over phi_i

            do i = 1, no_dofs_per_variable2(ieqn)

              face_residual_m(ieqn, i) = face_residual_m(ieqn, i) + integral_weighting(qk)*( &
                                         Re*(nflxsoln(ieqn, qk)*phi2(ieqn, qk, i)) &
               + (-0.5_db*1.000_db*(cal_bdry_sym_stress_u(gradient_uh1(:, qk, :), face_normals(:, qk), ieqn, problem_dim, no_pdes) &
                                  + cal_bdry_sym_stress_u(gradient_uh2(:, qk, :), face_normals(:, qk), ieqn, problem_dim, no_pdes) &
                                                              ) &
                                            + full_dispenal*(interpolant_uh1(ieqn, qk) - interpolant_uh2(ieqn, qk))) &
                                         *phi2(ieqn, qk, i) &
                                        + 0.5_db*1.000_db*cal_bdry_sym_stress_v(grad_phi2(ieqn, qk, :, i), interpolant_uh1(:, qk), &
                                                          interpolant_uh2(:, qk), face_normals(:, qk), ieqn, problem_dim, no_pdes) &
                                         + 0.5_db*(interpolant_uh1(problem_dim + 1, qk) + interpolant_uh2(problem_dim + 1, qk))* &
                                         face_normals(ieqn, qk)*phi2(ieqn, qk, i) &
                                         )

            end do

          end do

        end do

! =======================================================================================
! Divergence free constraint

        do qk = 1, no_quad_points

! Loop over phi_i

          do i = 1, no_dofs_per_variable1(problem_dim + 1)

            face_residual_p(no_pdes, i) = face_residual_p(no_pdes, i) + integral_weighting(qk)* &
                                          0.5_db*dot_product(interpolant_uh1(1:problem_dim, qk) &
                                                             - interpolant_uh2(1:problem_dim, qk), face_normals(:, qk)) &
                                          *phi1(problem_dim + 1, qk, i)

          end do

! Loop over phi_i

          do i = 1, no_dofs_per_variable2(problem_dim + 1)

            face_residual_m(no_pdes, i) = face_residual_m(no_pdes, i) + integral_weighting(qk)* &
                                          0.5_db*dot_product(interpolant_uh1(1:problem_dim, qk) &
                                                             - interpolant_uh2(1:problem_dim, qk), face_normals(:, qk)) &
                                          *phi2(problem_dim + 1, qk, i)

          end do

        end do

! =======================================================================================

      end if

    end associate

  end subroutine ns_dg_face_nonlin_residual

!--------------------------------------------------------------------
  subroutine ns_jac_element_mat(element_matrix, &
                                mesh_data, soln_data, facet_data, fe_basis_info)
!--------------------------------------------------------------------
    use problem_options
    use coupled_data_storage

    include 'assemble_jac_matrix_element.h'

! Local variables

    logical :: bdf_semi_implicit
    integer :: qk, i, j, ieqn, ivar
    real(db) :: diff_terms, gradgradterm, pgradv, qgradu, convection_term, brinkman_term
    real(db), dimension(facet_data%no_quad_points) :: eff_permeability_reciprocal
    real(db), dimension(facet_data%problem_dim, facet_data%problem_dim, &
                        facet_data%problem_dim, facet_data%no_quad_points) :: fluxes_prime
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points) :: interpolant_uh
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points, &
                        facet_data%problem_dim) :: gradient_uh
    real(db), dimension(facet_data%dim_soln_coeff, facet_data%no_quad_points, &
                        facet_data%problem_dim, maxval(facet_data%no_dofs_per_variable)) :: grad_phi
    real(db), dimension(facet_data%dim_soln_coeff, facet_data%no_quad_points, &
                        maxval(facet_data%no_dofs_per_variable)) :: phi
    real(db) :: div_u

    integer :: prev_sol, bdf_order
    real(db) :: bdf_scaling_factor
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points, &
                        facet_data%scheme_user_data%bdf_order) :: uh_previous_time_step
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points) :: &
      soln_extrap

    associate ( &
      dim_soln_coeff => facet_data%dim_soln_coeff, &
      no_pdes => facet_data%no_pdes, &
      problem_dim => facet_data%problem_dim, &
      no_quad_points => facet_data%no_quad_points, &
      global_points_ele => facet_data%global_points, &
      integral_weighting => facet_data%integral_weighting, &
      element_number => facet_data%element_number, &
      element_region_id => facet_data%element_region_id, &
      no_dofs_per_variable => facet_data%no_dofs_per_variable, &
      global_dof_numbers => facet_data%global_dof_numbers, &
      scheme_user_data => facet_data%scheme_user_data)

      bdf_order = scheme_user_data%bdf_order
      bdf_semi_implicit = scheme_user_data%bdf_semi_implicit
      bdf_scaling_factor = scheme_user_data%bdf_scaling_factor

      select type (scheme_user_data)
      type is (coupled_user_data)
        do prev_sol = 1, bdf_order
          call compute_uh_with_basis_fns_pts_from_array(uh_previous_time_step(:, :, prev_sol), &
                                                        no_pdes, no_quad_points, dim_soln_coeff, no_dofs_per_variable, &
                                                        global_dof_numbers, fe_basis_info%basis_element, soln_data, &
                                                        soln_data%no_dofs, fe_data_struct(1)%soln_data_prev(prev_sol)%soln_values)
        end do
      class default
        do prev_sol = 1, bdf_order
          call compute_uh_with_basis_fns_pts_from_array(uh_previous_time_step(:, :, prev_sol), &
                                                        no_pdes, no_quad_points, dim_soln_coeff, no_dofs_per_variable, &
                                                        global_dof_numbers, fe_basis_info%basis_element, soln_data, &
                                                        soln_data%no_dofs, scheme_user_data%uh_previous_time_step_ms(:, prev_sol))
        end do
      end select

      element_matrix = 0.0_db

      do qk = 1, no_quad_points
        interpolant_uh(:, qk) = uh_element(fe_basis_info, no_pdes, qk)

        do i = 1, no_pdes
          gradient_uh(i, qk, 1:problem_dim) = grad_uh_element(fe_basis_info, problem_dim, i, qk, 1)
        end do

        call jac_adv_flux_bdf(fluxes_prime(:, :, :, qk), interpolant_uh(:, qk), problem_dim, no_pdes, &
                              bdf_order, bdf_semi_implicit, uh_previous_time_step(:, qk, :))

        call eval_tanh_smooth(problem_dim, element_number, global_points_ele(:, qk), &
                              Da_reciprocal, eff_permeability_reciprocal(qk))

      end do

! Calculate Basis Functions

      do i = 1, dim_soln_coeff
        grad_phi(i, 1:no_quad_points, 1:problem_dim, 1:no_dofs_per_variable(i)) = fe_basis_info%basis_element &
                                        %deriv_basis_fns(i)%grad_data(1:no_quad_points, 1:problem_dim, 1:no_dofs_per_variable(i), 1)
        phi(i, 1:no_quad_points, 1:no_dofs_per_variable(i)) = fe_basis_info%basis_element%basis_fns(i) &
                                                              %fem_basis_fns(1:no_quad_points, 1:no_dofs_per_variable(i), 1)
      end do

! Loop over quadrature points

      do qk = 1, no_quad_points

        div_u = 0.0_db

        do ieqn = 1, problem_dim

          div_u = div_u + gradient_uh(ieqn, qk, ieqn)

        end do

! Loop over the equations

        do ieqn = 1, no_pdes

! Loop over phi_i

          do i = 1, no_dofs_per_variable(ieqn)

! Loop over the variables

            do ivar = 1, no_pdes

! Loop over phi_j

              do j = 1, no_dofs_per_variable(ivar)

                convection_term = 0.0_db
                brinkman_term = 0.0_db
                if (ieqn <= problem_dim .and. ivar <= problem_dim) then

                  convection_term = dot_product( &
                                    fluxes_prime(1:problem_dim, ieqn, ivar, qk), &
                                    grad_phi(ieqn, qk, 1:problem_dim, i))*phi(ivar, qk, j)

                  if (ivar == ieqn) then
                    !convection_term = convection_term+(div_u*phi(ivar,qk,j)+&
                    !     grad_phi(ivar,qk,ivar,j)*interpolant_uh(ieqn,qk))*&
                    !     phi(ieqn,qk,i)/2.0_db

                    ! Time term
                    ! Apto takes -ve Jac. so sign should be same as for F here
                    ! partial u_i / partial u_j = krocknecker_delta_i^j
                    convection_term = convection_term &
                                      - bdf_scaling_factor*phi(ieqn, qk, i)*phi(ivar, qk, j)

                    brinkman_term = brinkman_term &
                                    - phi(ieqn, qk, i)*phi(ivar, qk, j)

                    !else
                    !convection_term = convection_term+grad_phi(ivar,qk,ivar,j)*&
                    !     interpolant_uh(ieqn,qk)*&
                    !     phi(ieqn,qk,i)/2.0_db
                  end if

                end if

                gradgradterm = 0.5_db*cal_sym_stress_deriv(grad_phi(:, qk, :, j), &
                                                           grad_phi(:, qk, :, i), ivar, ieqn, problem_dim, no_pdes)

                pgradv = cal_gradterm(phi(ivar, qk, j), &
                                      grad_phi(ieqn, qk, :, i), ivar, ieqn, problem_dim, no_pdes)

                qgradu = cal_gradterm(phi(ieqn, qk, i), &
                                      grad_phi(ivar, qk, :, j), ieqn, ivar, problem_dim, no_pdes)

                diff_terms = 1.000_db*gradgradterm - pgradv + qgradu

                convection_term = Re*convection_term

                brinkman_term = 1.000_db*eff_permeability_reciprocal(qk)*brinkman_term

                element_matrix(ieqn, ivar, i, j) = element_matrix(ieqn, ivar, i, j) &
                                                   + integral_weighting(qk)*(diff_terms - convection_term - brinkman_term)

              end do

            end do

          end do

        end do

      end do

    end associate

  end subroutine ns_jac_element_mat

!--------------------------------------------------------------------
  subroutine ns_dg_jac_face_mat(face_matrix_pp, face_matrix_pm, &
                                face_matrix_mp, face_matrix_mm, mesh_data, soln_data, &
                                facet_data, fe_basis_info)
!--------------------------------------------------------------------

    use problem_options
    use coupled_data_storage
    use aptofem_fe_matrix_assembly

    include 'assemble_jac_matrix_int_bdry_face.h'

! Local variables

    integer :: qk, i, j, ieqn, ivar
    real(db) :: full_dispenal, diff_terms, pvterm, quterm, &
                convective_flux_1, convective_flux_2, convection_terms, robin_bc_term, bf_stab
    real(db), dimension(facet_data%problem_dim, facet_data%problem_dim, &
                        facet_data%problem_dim, facet_data%no_quad_points) :: &
      fluxes_prime1, fluxes_prime2
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points) :: uloc
    real(db), dimension(facet_data%no_quad_points) :: alpha
    real(db), dimension(facet_data%problem_dim, facet_data%problem_dim, &
                        facet_data%no_quad_points) :: boundary_jacobian
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points) :: &
      interpolant_uh1, interpolant_uh2
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points, &
                        facet_data%problem_dim) :: grad_uh1, grad_uh2
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points, &
                        facet_data%problem_dim, maxval(facet_data%no_dofs_per_variable1)) :: grad_phi1
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points, &
                        facet_data%problem_dim, maxval(facet_data%no_dofs_per_variable2)) :: grad_phi2
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points, &
                        maxval(facet_data%no_dofs_per_variable1)) :: phi1
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points, &
                        maxval(facet_data%no_dofs_per_variable2)) :: phi2
    real(db), dimension(facet_data%dim_soln_coeff) :: dispenal_local

    logical :: bdf_semi_implicit
    integer :: prev_sol, bdf_order
    real(db) :: bdf_scaling_factor, current_time, tang_pen, udotn_factor_bdf
    real(db), dimension(facet_data%no_quad_points) :: bf_stab_term
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points, &
                        facet_data%scheme_user_data%bdf_order) :: uh_previous_time_step1
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points, &
                        facet_data%scheme_user_data%bdf_order) :: uh_previous_time_step2
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points, &
                        facet_data%problem_dim, facet_data%scheme_user_data%bdf_order) :: &
      gradient_uh_previous_time_step1
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points, &
                        facet_data%problem_dim, facet_data%scheme_user_data%bdf_order) :: &
      gradient_uh_previous_time_step2
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points) :: &
      soln_extrap1
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points) :: &
      soln_extrap2

    associate ( &
      dim_soln_coeff => facet_data%dim_soln_coeff, &
      no_pdes => facet_data%no_pdes, &
      problem_dim => facet_data%problem_dim, &
      no_quad_points => facet_data%no_quad_points, &
      global_points_face => facet_data%global_points, &
      integral_weighting => facet_data%integral_weighting, &
      face_number => facet_data%face_number, &
      neighbours => facet_data%neighbours, &
      interior_face_boundary_no => facet_data%interior_face_boundary_no, &
      face_element_region_ids => facet_data%face_element_region_ids, &
      bdry_face => facet_data%bdry_no, &
      no_dofs_per_variable1 => facet_data%no_dofs_per_variable1, &
      no_dofs_per_variable2 => facet_data%no_dofs_per_variable2, &
      face_normals => facet_data%face_normals, &
      dispenal => facet_data%dispenal, &
      scheme_user_data => facet_data%scheme_user_data)

      bdf_order = scheme_user_data%bdf_order
      bdf_semi_implicit = scheme_user_data%bdf_semi_implicit
      bdf_scaling_factor = scheme_user_data%bdf_scaling_factor
      current_time = scheme_user_data%current_time

      face_matrix_pp = 0.0_db
      face_matrix_pm = 0.0_db
      face_matrix_mp = 0.0_db
      face_matrix_mm = 0.0_db

      dispenal_local = dispenal

      full_dispenal = interior_penalty_parameter*dispenal_local(1)

      if (bdry_face > 0) then

        select type (scheme_user_data)
        type is (coupled_user_data)
          do prev_sol = 1, bdf_order
            call compute_uh_with_basis_fns_pts_from_array(uh_previous_time_step1(:, :, prev_sol), &
                                                        no_pdes, no_quad_points, dim_soln_coeff, facet_data%no_dofs_per_variable1, &
                                                          facet_data%global_dof_numbers1, fe_basis_info%basis_face1, soln_data, &
                                                          soln_data%no_dofs, fe_data_struct(1)%soln_data_prev(prev_sol)%soln_values)
          end do
        class default
          if (bdry_face > 300) then
            write (io_err, *) 'ERROR: bdry_face > 300 can only be used with coupled_user_data type'
            stop
          end if
          do prev_sol = 1, bdf_order
            call compute_uh_with_basis_fns_pts_from_array(uh_previous_time_step1(:, :, prev_sol), &
                                                        no_pdes, no_quad_points, dim_soln_coeff, facet_data%no_dofs_per_variable1, &
                                                          facet_data%global_dof_numbers1, fe_basis_info%basis_face1, soln_data, &
                                                          soln_data%no_dofs, scheme_user_data%uh_previous_time_step_ms(:, prev_sol))
          end do
        end select

! ====================================================================
! Boundary Face

! Calculate value of analytical solution at quadrature points

        do qk = 1, no_quad_points
          interpolant_uh1(:, qk) = uh_face1(fe_basis_info, no_pdes, qk)
          do i = 1, no_pdes
            grad_uh1(i, qk, 1:problem_dim) = grad_uh_face1(fe_basis_info, problem_dim, i, qk, 1)
          end do

          call compute_boundary_condition_bdf(interpolant_uh2(:, qk), &
                                              uloc(:, qk), interpolant_uh1(:, qk), face_normals(:, qk), &
                                              problem_dim, no_pdes, abs(bdry_face), global_points_face(:, qk), &
                                              scheme_user_data%current_time)
          call lax_friedrichs_jac_bdf_bdry(boundary_jacobian(:, :, qk), &
                                           fluxes_prime1(:, :, :, qk), fluxes_prime2(:, :, :, qk), alpha(qk), &
                             interpolant_uh1(:, qk), interpolant_uh2(:, qk), face_normals(:, qk), problem_dim, no_pdes, bdry_face, &
                                           bdf_order, bdf_semi_implicit, uh_previous_time_step1(:, qk, :))

          call bdf_extrapolate(soln_extrap1(:, qk), uh_previous_time_step1(:, qk, :), &
                               no_pdes, bdf_order)
          call normal_vel_stab(bf_stab_term(qk), &
                               soln_extrap1(:, qk), face_normals(:, qk), &
                               problem_dim, no_pdes, backflow_stab_parameter)
        end do

        do i = 1, no_pdes
          grad_phi1(i, 1:no_quad_points, 1:problem_dim, 1:no_dofs_per_variable1(i)) = fe_basis_info%basis_face1 &
                                                   %deriv_basis_fns(i)%grad_data(1:no_quad_points, :, 1:no_dofs_per_variable1(i), 1)
          phi1(i, 1:no_quad_points, 1:no_dofs_per_variable1(i)) = fe_basis_info%basis_face1%basis_fns(i) &
                                                                  %fem_basis_fns(1:no_quad_points, 1:no_dofs_per_variable1(i), 1)
        end do

        if (abs(bdry_face) <= 100) then ! Dirichlet boundary

! Loop over quadrature points

          do qk = 1, no_quad_points

! Loop over the equations

            do ieqn = 1, no_pdes

! Loop over phi_i

              do i = 1, no_dofs_per_variable1(ieqn)

! Loop over the variables

                do ivar = 1, no_pdes

! Loop over phi_j

                  do j = 1, no_dofs_per_variable1(ivar)

                    convective_flux_1 = 0.0_db
                    convective_flux_2 = 0.0_db

                    if (ieqn <= problem_dim .and. ivar <= problem_dim) then
                      convective_flux_1 = 0.5_db*dot_product( &
                                          fluxes_prime1(1:problem_dim, ieqn, ivar, qk), &
                                          face_normals(1:problem_dim, qk))
                      if (ieqn == ivar) then
                        convective_flux_1 = convective_flux_1 + 0.5_db*alpha(qk)
                      end if
                      convective_flux_2 = boundary_jacobian(ieqn, ivar, qk)

                    end if

                    convection_terms = (convective_flux_1 + convective_flux_2) &
                                       *phi1(ivar, qk, j)*phi1(ieqn, qk, i)

                    diff_terms = cal_stress_terms_bdry_face(grad_phi1(ivar, qk, :, j), &
                                                            grad_phi1(ieqn, qk, :, i), phi1(ivar, qk, j), &
                                                            phi1(ieqn, qk, i), face_normals(:, qk), full_dispenal, &
                                                            ivar, ieqn, problem_dim, no_pdes)

                    pvterm = cal_grad_terms_bdry_face(phi1(ivar, qk, j), &
                                                      phi1(ieqn, qk, i), face_normals(:, qk), ivar, ieqn, &
                                                      problem_dim, no_pdes)

                    quterm = cal_grad_terms_bdry_face(phi1(ieqn, qk, i), &
                                                      phi1(ivar, qk, j), face_normals(:, qk), ieqn, ivar, &
                                                      problem_dim, no_pdes)

                    face_matrix_pp(ieqn, ivar, i, j) = face_matrix_pp(ieqn, ivar, i, j) &
                                                       + integral_weighting(qk)*(1.000_db*diff_terms + pvterm - quterm &
                                                                                 + Re*convection_terms)

                  end do

                end do

              end do

            end do

          end do

        elseif (abs(bdry_face) <= 200) then ! Neumann Boundary

! Loop over quadrature points

          do qk = 1, no_quad_points

! Loop over the equations

            do ieqn = 1, no_pdes

! Loop over phi_i

              do i = 1, no_dofs_per_variable1(ieqn)

! Loop over the variables

                do ivar = 1, no_pdes

! Loop over phi_j

                  do j = 1, no_dofs_per_variable1(ivar)

                    convective_flux_1 = 0.0_db
                    convective_flux_2 = 0.0_db
                    bf_stab = 0.0_db

                    if (ieqn <= problem_dim .and. ivar <= problem_dim) then
                      convective_flux_1 = 0.5_db*dot_product( &
                                          fluxes_prime1(1:problem_dim, ieqn, ivar, qk), &
                                          face_normals(1:problem_dim, qk))
                      if (ieqn == ivar) then
                        convective_flux_1 = convective_flux_1 + 0.5_db*alpha(qk)
                        bf_stab = bf_stab + bf_stab_term(qk)
                      end if
                      convective_flux_2 = boundary_jacobian(ieqn, ivar, qk)

                    end if

                    convection_terms = (convective_flux_1 + convective_flux_2)

                    face_matrix_pp(ieqn, ivar, i, j) = face_matrix_pp(ieqn, ivar, i, j) &
                                                       + integral_weighting(qk) &
                                                       *Re*(bf_stab*phi1(ivar, qk, j)*phi1(ieqn, qk, i) &
                                                            + convection_terms*phi1(ivar, qk, j)*phi1(ieqn, qk, i))

                  end do

                end do

              end do

            end do

          end do

        end if

! End of boundary face
! ====================================================================

      else

! coupled_user_data type implies fsi_bdf_timestepping, otherwise bdf_timestepping
        select type (scheme_user_data)
        type is (coupled_user_data)
          do prev_sol = 1, bdf_order
            call compute_uh_with_basis_fns_pts_from_array(uh_previous_time_step1(:, :, prev_sol), &
                                                        no_pdes, no_quad_points, dim_soln_coeff, facet_data%no_dofs_per_variable1, &
                                                          facet_data%global_dof_numbers1, fe_basis_info%basis_face1, soln_data, &
                                                          soln_data%no_dofs, fe_data_struct(1)%soln_data_prev(prev_sol)%soln_values)
            call compute_uh_with_basis_fns_pts_from_array(uh_previous_time_step2(:, :, prev_sol), &
                                                        no_pdes, no_quad_points, dim_soln_coeff, facet_data%no_dofs_per_variable2, &
                                                          facet_data%global_dof_numbers2, fe_basis_info%basis_face2, soln_data, &
                                                          soln_data%no_dofs, fe_data_struct(1)%soln_data_prev(prev_sol)%soln_values)
          end do
        class default
          do prev_sol = 1, bdf_order
            call compute_uh_with_basis_fns_pts_from_array(uh_previous_time_step1(:, :, prev_sol), &
                                                        no_pdes, no_quad_points, dim_soln_coeff, facet_data%no_dofs_per_variable1, &
                                                          facet_data%global_dof_numbers1, fe_basis_info%basis_face1, soln_data, &
                                                          soln_data%no_dofs, scheme_user_data%uh_previous_time_step_ms(:, prev_sol))
            call compute_uh_with_basis_fns_pts_from_array(uh_previous_time_step2(:, :, prev_sol), &
                                                        no_pdes, no_quad_points, dim_soln_coeff, facet_data%no_dofs_per_variable2, &
                                                          facet_data%global_dof_numbers2, fe_basis_info%basis_face2, soln_data, &
                                                          soln_data%no_dofs, scheme_user_data%uh_previous_time_step_ms(:, prev_sol))
          end do
        end select

! ====================================================================
! Interior Face

        do qk = 1, no_quad_points
          interpolant_uh1(:, qk) = uh_face1(fe_basis_info, no_pdes, qk)
          interpolant_uh2(:, qk) = uh_face2(fe_basis_info, no_pdes, qk)
          do i = 1, no_pdes
            grad_uh1(i, qk, 1:problem_dim) = grad_uh_face1(fe_basis_info, problem_dim, i, qk, 1)
            grad_uh2(i, qk, 1:problem_dim) = grad_uh_face2(fe_basis_info, problem_dim, i, qk, 1)
          end do

          call bdf_extrapolate(soln_extrap1(:, qk), uh_previous_time_step1(:, qk, :), &
                               no_pdes, bdf_order)
          call bdf_extrapolate(soln_extrap2(:, qk), uh_previous_time_step2(:, qk, :), &
                               no_pdes, bdf_order)

        end do

        do qk = 1, no_quad_points
          alpha(qk) = 2.0_db*max(abs(dot_product(soln_extrap1(1:problem_dim, qk), face_normals(:, qk))), &
                                 abs(dot_product(soln_extrap2(1:problem_dim, qk), face_normals(:, qk))))
          call jac_adv_flux_extrap(fluxes_prime1(:, :, :, qk), problem_dim, no_pdes, &
                                   soln_extrap1(:, qk))
          call jac_adv_flux_extrap(fluxes_prime2(:, :, :, qk), problem_dim, no_pdes, &
                                   soln_extrap2(:, qk))
        end do

        do i = 1, no_pdes
          grad_phi1(i, 1:no_quad_points, 1:problem_dim, 1:no_dofs_per_variable1(i)) = fe_basis_info%basis_face1 &
                                                   %deriv_basis_fns(i)%grad_data(1:no_quad_points, :, 1:no_dofs_per_variable1(i), 1)
          grad_phi2(i, 1:no_quad_points, 1:problem_dim, 1:no_dofs_per_variable2(i)) = fe_basis_info%basis_face2 &
                                                   %deriv_basis_fns(i)%grad_data(1:no_quad_points, :, 1:no_dofs_per_variable2(i), 1)
          phi1(i, 1:no_quad_points, 1:no_dofs_per_variable1(i)) = fe_basis_info%basis_face1%basis_fns(i) &
                                                                  %fem_basis_fns(1:no_quad_points, 1:no_dofs_per_variable1(i), 1)
          phi2(i, 1:no_quad_points, 1:no_dofs_per_variable2(i)) = fe_basis_info%basis_face2%basis_fns(i) &
                                                                  %fem_basis_fns(1:no_quad_points, 1:no_dofs_per_variable2(i), 1)
        end do

        do qk = 1, no_quad_points

          do ieqn = 1, no_pdes

            do ivar = 1, no_pdes

              convective_flux_1 = 0.0_db
              convective_flux_2 = 0.0_db

              if (ieqn <= problem_dim .and. ivar <= problem_dim) then
                convective_flux_1 = 0.5_db*dot_product( &
                                    fluxes_prime1(1:problem_dim, ieqn, ivar, qk), &
                                    face_normals(1:problem_dim, qk))
                convective_flux_2 = 0.5_db*dot_product( &
                                    fluxes_prime2(1:problem_dim, ieqn, ivar, qk), &
                                    face_normals(1:problem_dim, qk))
                if (ieqn == ivar) then
                  convective_flux_1 = convective_flux_1 + 0.5_db*alpha(qk)
                  convective_flux_2 = convective_flux_2 - 0.5_db*alpha(qk)
                end if

              end if

! Loop over phi_i

              do i = 1, no_dofs_per_variable1(ieqn)

! u^+ v^+ w.r.t. ele1

! Loop over phi_j

                do j = 1, no_dofs_per_variable1(ivar)

                  diff_terms = cal_stress_terms_int_face(grad_phi1(ivar, qk, :, j), &
                                                         grad_phi1(ieqn, qk, :, i), phi1(ivar, qk, j), &
                                                         phi1(ieqn, qk, i), face_normals(:, qk), full_dispenal, &
                                                         ivar, ieqn, problem_dim, no_pdes)

                  pvterm = 0.5_db*cal_grad_terms_bdry_face(phi1(ivar, qk, j), &
                                                           phi1(ieqn, qk, i), face_normals(:, qk), ivar, ieqn, &
                                                           problem_dim, no_pdes)

                  quterm = 0.5_db*cal_grad_terms_bdry_face(phi1(ieqn, qk, i), &
                                                           phi1(ivar, qk, j), face_normals(:, qk), ieqn, ivar, &
                                                           problem_dim, no_pdes)

                  face_matrix_pp(ieqn, ivar, i, j) = face_matrix_pp(ieqn, ivar, i, j) &
                                                     + integral_weighting(qk)*(1.000_db*diff_terms + pvterm - quterm &
                                                                         + Re*convective_flux_1*phi1(ieqn, qk, i)*phi1(ivar, qk, j))

                end do

! u^- v^+ w.r.t. ele1

! Loop over phi_j

                do j = 1, no_dofs_per_variable2(ivar)

                  diff_terms = cal_stress_terms_int_face(grad_phi2(ivar, qk, :, j), &
                                                         grad_phi1(ieqn, qk, :, i), -phi2(ivar, qk, j), &
                                                         phi1(ieqn, qk, i), face_normals(:, qk), full_dispenal, &
                                                         ivar, ieqn, problem_dim, no_pdes)

                  pvterm = 0.5_db*cal_grad_terms_bdry_face(phi2(ivar, qk, j), &
                                                           phi1(ieqn, qk, i), face_normals(:, qk), ivar, ieqn, &
                                                           problem_dim, no_pdes)

                  quterm = -0.5_db*cal_grad_terms_bdry_face(phi1(ieqn, qk, i), &
                                                            phi2(ivar, qk, j), face_normals(:, qk), ieqn, ivar, &
                                                            problem_dim, no_pdes)

                  face_matrix_mp(ieqn, ivar, i, j) = face_matrix_mp(ieqn, ivar, i, j) &
                                                     + integral_weighting(qk)*(1.000_db*diff_terms + pvterm - quterm &
                                                                         + Re*convective_flux_2*phi1(ieqn, qk, i)*phi2(ivar, qk, j))

                end do

              end do

! Loop over phi_i

              do i = 1, no_dofs_per_variable2(ieqn)

! u^+ v^- w.r.t. ele1

! Loop over phi_j

                do j = 1, no_dofs_per_variable1(ivar)

                  diff_terms = cal_stress_terms_int_face(grad_phi1(ivar, qk, :, j), &
                                                         grad_phi2(ieqn, qk, :, i), phi1(ivar, qk, j), &
                                                         -phi2(ieqn, qk, i), face_normals(:, qk), full_dispenal, &
                                                         ivar, ieqn, problem_dim, no_pdes)

                  pvterm = -0.5_db*cal_grad_terms_bdry_face(phi1(ivar, qk, j), &
                                                            phi2(ieqn, qk, i), face_normals(:, qk), ivar, ieqn, &
                                                            problem_dim, no_pdes)

                  quterm = 0.5_db*cal_grad_terms_bdry_face(phi2(ieqn, qk, i), &
                                                           phi1(ivar, qk, j), face_normals(:, qk), ieqn, ivar, &
                                                           problem_dim, no_pdes)

                  face_matrix_pm(ieqn, ivar, i, j) = face_matrix_pm(ieqn, ivar, i, j) &
                                                     + integral_weighting(qk)*(1.000_db*diff_terms + pvterm - quterm &
                                                                         - Re*convective_flux_1*phi2(ieqn, qk, i)*phi1(ivar, qk, j))

                end do

! u^- v^- w.r.t. ele1

! Loop over phi_j

                do j = 1, no_dofs_per_variable2(ivar)

                  diff_terms = cal_stress_terms_int_face(grad_phi2(ivar, qk, :, j), &
                                                         grad_phi2(ieqn, qk, :, i), -phi2(ivar, qk, j), &
                                                         -phi2(ieqn, qk, i), face_normals(:, qk), full_dispenal, &
                                                         ivar, ieqn, problem_dim, no_pdes)

                  pvterm = -0.5_db*cal_grad_terms_bdry_face(phi2(ivar, qk, j), &
                                                            phi2(ieqn, qk, i), face_normals(:, qk), ivar, ieqn, &
                                                            problem_dim, no_pdes)

                  quterm = -0.5_db*cal_grad_terms_bdry_face(phi2(ieqn, qk, i), &
                                                            phi2(ivar, qk, j), face_normals(:, qk), ieqn, ivar, &
                                                            problem_dim, no_pdes)

                  face_matrix_mm(ieqn, ivar, i, j) = face_matrix_mm(ieqn, ivar, i, j) &
                                                     + integral_weighting(qk)*(1.000_db*diff_terms + pvterm - quterm &
                                                                         - Re*convective_flux_2*phi2(ieqn, qk, i)*phi2(ivar, qk, j))

                end do

              end do

            end do

          end do

        end do

! End of interior face
! ====================================================================
        
      end if

    end associate

  end subroutine ns_dg_jac_face_mat

end module navier_stokes_jac_matrix_residual
