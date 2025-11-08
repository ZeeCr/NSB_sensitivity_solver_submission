module adr_av_jac_matrix_residual

  use param

  use coupled_data_storage
  use problem_options
  use geometry_routines

  real(db), parameter :: C_AV = 1.0_db/1.0_db
  real(db), parameter :: B_AV = 2.0_db/10.0_db

contains

!--------------------------------------------------------------------
!> Residual for ADR
!--------------------------------------------------------------------
  subroutine adr_av_element_nonlin_residual(element_rhs, &
                                            mesh_data, soln_data, facet_data, fe_basis_info)
!--------------------------------------------------------------------

    include 'assemble_residual_element.h'

! Local variables

    integer :: qk, i, j, ieqn
    real(db), dimension(facet_data%no_quad_points) :: eff_Dm
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points) &
      :: floc
    real(db) :: div_u
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points) :: &
      interpolant_uh
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points, &
                        facet_data%problem_dim) :: gradient_uh
    real(db), dimension(facet_data%dim_soln_coeff, facet_data%no_quad_points, &
                        facet_data%problem_dim, maxval(facet_data%no_dofs_per_variable)) :: grad_phi
    real(db), dimension(facet_data%dim_soln_coeff, facet_data%no_quad_points, &
                        maxval(facet_data%no_dofs_per_variable)) :: phi

    integer :: dim_soln_coeff_fluid, no_pdes_fluid, problem_dim_fluid
    real(db) :: adv_terms, diff_terms, reaction_terms, residual, residual_terms, &
                ele_diam, coeff_AV
    real(db), dimension(:, :), allocatable :: adv_vel
    real(db), dimension(:, :, :), allocatable :: &
      grad_adv_vel

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

      associate ( &
        soln_data_fluid => fe_data_struct(1)%external_soln, &
        mesh_data_fluid => fe_data_struct(1)%external_mesh)
        dim_soln_coeff_fluid = get_dim_soln_coeff(soln_data_fluid)
        no_pdes_fluid = get_no_pdes(soln_data_fluid)
        problem_dim_fluid = get_problem_dim(mesh_data_fluid)
        allocate (adv_vel(no_pdes_fluid, facet_data%no_quad_points), &
                  grad_adv_vel(no_pdes_fluid, problem_dim_fluid, facet_data%no_quad_points))

        element_rhs = 0.0_db

        div_u = 0.0_db
        residual = 0.0_db
        ele_diam = compute_h_iso(mesh_data, element_number)
        coeff_AV = C_AV*(ele_diam**(2.0_db - B_AV))

        do qk = 1, no_quad_points

          call compute_uh_gradient_uh_glob_pt(adv_vel(:, qk), grad_adv_vel(:, :, qk), &
                                              no_pdes_fluid, element_number, &
                                              global_points_ele(:, qk), problem_dim_fluid, mesh_data_fluid, soln_data_fluid)
          adv_vel(problem_dim_fluid + 1, qk) = 0.0_db
          grad_adv_vel(problem_dim_fluid + 1, :, qk) = 0.0_db

          interpolant_uh(1:no_pdes, qk) = uh_element(fe_basis_info, no_pdes, qk)
          do i = 1, no_pdes
            gradient_uh(i, qk, 1:problem_dim) = grad_uh_element(fe_basis_info, problem_dim, i, qk, 1)
          end do
          call forcing_fn(floc(:, qk), global_points_ele(:, qk), problem_dim, no_pdes, 0.0_db)

          call eval_tanh_smooth(problem_dim, element_number, global_points_ele(:, qk), &
                                Dm, eff_Dm(qk))

        end do

! Calculate Basis Functions

        do i = 1, dim_soln_coeff
          grad_phi(i, 1:no_quad_points, 1:problem_dim, 1:no_dofs_per_variable(i)) = fe_basis_info%basis_element &
                                        %deriv_basis_fns(i)%grad_data(1:no_quad_points, 1:problem_dim, 1:no_dofs_per_variable(i), 1)
          phi(i, 1:no_quad_points, 1:no_dofs_per_variable(i)) = fe_basis_info%basis_element%basis_fns(i) &
                                                                %fem_basis_fns(1:no_quad_points, 1:no_dofs_per_variable(i), 1)
        end do

! Momentum Equations

        do ieqn = 1, no_pdes

! Loop over quadrature points

          do qk = 1, no_quad_points

            do j = 1, problem_dim_fluid
              div_u = div_u + grad_adv_vel(j, j, qk)
            end do

! Loop over phi_i

            do i = 1, no_dofs_per_variable(ieqn)

              residual = &
                abs( &
                interpolant_uh(ieqn, qk)*div_u &
                + dot_product(adv_vel(1:problem_dim_fluid, qk), gradient_uh(ieqn, qk, :)) &
                + eff_Dm(qk)*interpolant_uh(ieqn, qk) &
                )

              adv_terms = interpolant_uh(ieqn, qk)*dot_product( &
                          adv_vel(1:problem_dim_fluid, qk), grad_phi(ieqn, qk, :, i))

              diff_terms = (1.0_db/Pe)*( &
                           -dot_product( &
                           gradient_uh(ieqn, qk, :), grad_phi(ieqn, qk, :, i)) &
                           )

              reaction_terms = -eff_Dm(qk)*interpolant_uh(ieqn, qk)*phi(ieqn, qk, i)

              residual_terms = &
                -coeff_AV* &
                residual*dot_product( &
                gradient_uh(ieqn, qk, :), grad_phi(ieqn, qk, :, i) &
                )

              element_rhs(ieqn, i) = element_rhs(ieqn, i) &
                                     + integral_weighting(qk)*( &
                                     adv_terms + diff_terms + reaction_terms &
                                     + residual_terms + floc(ieqn, qk)*phi(ieqn, qk, i) &
                                     )

            end do
          end do
        end do

        deallocate (adv_vel, grad_adv_vel)

      end associate

    end associate

  end subroutine adr_av_element_nonlin_residual

!--------------------------------------------------------------------
!> Face residual for ADR
!  --------------------------------------------------------------
  subroutine adr_av_dg_face_nonlin_residual(face_residual_p, face_residual_m, &
                                            mesh_data, soln_data, facet_data, fe_basis_info)
!--------------------------------------------------------------------
    use aptofem_fe_matrix_assembly

    include 'assemble_residual_int_bdry_face.h'

! Local variables

    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points) &
      :: uloc
    real(db), dimension(facet_data%no_pdes, &
                        facet_data%no_quad_points) :: unloc
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

    integer :: dim_soln_coeff_fluid, no_pdes_fluid, problem_dim_fluid
    logical, dimension(facet_data%no_quad_points) :: flow_out
    real(db) :: adv_terms, diff_terms, upwind_val_int_face, tol, gR, robin_term
    real(db), dimension(facet_data%no_quad_points) :: normal_vel
    real(db), dimension(:, :), allocatable :: &
      adv_vel1, adv_vel2

    tol = 1.0d-6

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

      associate ( &
        soln_data_fluid => fe_data_struct(1)%external_soln, &
        mesh_data_fluid => fe_data_struct(1)%external_mesh)
        dim_soln_coeff_fluid = get_dim_soln_coeff(soln_data_fluid)
        no_pdes_fluid = get_no_pdes(soln_data_fluid)
        problem_dim_fluid = get_problem_dim(mesh_data_fluid)
        allocate (adv_vel1(no_pdes_fluid, facet_data%no_quad_points), &
                  adv_vel2(no_pdes_fluid, facet_data%no_quad_points))

        face_residual_p = 0.0_db
        face_residual_m = 0.0_db

        dispenal_local = dispenal

        full_dispenal = interior_penalty_parameter*dispenal_local(1)*Pe

        if (bdry_face > 0) then

          do qk = 1, no_quad_points
            interpolant_uh1(1:no_pdes, qk) = uh_face1(fe_basis_info, no_pdes, qk)
            do i = 1, no_pdes
              gradient_uh1(i, qk, 1:problem_dim) = grad_uh_face1(fe_basis_info, problem_dim, i, qk, 1)
            end do

            do i = 1, no_pdes
              grad_phi1(i, 1:no_quad_points, 1:problem_dim, 1:no_dofs_per_variable1(i)) = fe_basis_info%basis_face1 &
                                                   %deriv_basis_fns(i)%grad_data(1:no_quad_points, :, 1:no_dofs_per_variable1(i), 1)
              phi1(i, 1:no_quad_points, 1:no_dofs_per_variable1(i)) = fe_basis_info%basis_face1%basis_fns(i) &
                                                                     %fem_basis_fns(1:no_quad_points, 1:no_dofs_per_variable1(i), 1)
            end do

            call compute_uh_glob_pt(adv_vel1(:, qk), &
                                    no_pdes_fluid, neighbours(1), &
                                    global_points_face(:, qk), problem_dim_fluid, mesh_data_fluid, soln_data_fluid)
            adv_vel1(problem_dim_fluid + 1, qk) = 0.0_db

            normal_vel(qk) = dot_product(adv_vel1(1:problem_dim_fluid, qk), face_normals(:, qk))
            if (normal_vel(qk) > -tol) then
              flow_out(qk) = .TRUE.
            else
              flow_out(qk) = .FALSE.
            end if

          end do

          do ieqn = 1, no_pdes
            do qk = 1, no_quad_points

              if (flow_out(qk) .AND. &
                  bdry_face >= 101 .AND. bdry_face <= 200) then

! Neumann boundary

! =======================================================================================

                call neumann_bc(unloc(:, qk), global_points_face(:, qk), problem_dim, bdry_face, 0.0_db)

! Loop over phi_i

                do i = 1, no_dofs_per_variable1(ieqn)

                  adv_terms = -interpolant_uh1(ieqn, qk)*phi1(ieqn, qk, i)*normal_vel(qk)

                  diff_terms = (1.0_db/Pe)*( &
                               unloc(ieqn, qk)*phi1(ieqn, qk, i) &
                               )

                  ! RHS
                  face_residual_p(ieqn, i) = face_residual_p(ieqn, i) &
                                             + integral_weighting(qk)*( &
                                             adv_terms + diff_terms &
                                             )

                end do

                ! Inlets, i.e. Dirichlet boundary minus no-slip boundary
              else if (bdry_face >= 2 .AND. bdry_face <= 100) then

                call anal_soln(uloc(:, qk), global_points_face(:, qk), problem_dim, no_pdes, bdry_face, 0.0_db)

                ! Loop over phi_i

                do i = 1, no_dofs_per_variable1(ieqn)

                  adv_terms = -uloc(ieqn, qk)*phi1(ieqn, qk, i)*normal_vel(qk)

                  diff_terms = (1.0_db/Pe)*( &
                               -full_dispenal*(interpolant_uh1(ieqn, qk) - uloc(ieqn, qk))*phi1(ieqn, qk, i) &
                               + dot_product(gradient_uh1(ieqn, qk, :), face_normals(:, qk))*phi1(ieqn, qk, i) &
                               + dot_product(grad_phi1(ieqn, qk, :, i), face_normals(:, qk))* &
                               (interpolant_uh1(ieqn, qk) - uloc(ieqn, qk)) &
                               )

                  ! RHS
                  face_residual_p(ieqn, i) = face_residual_p(ieqn, i) &
                                             + integral_weighting(qk)*( &
                                             adv_terms + diff_terms)

                end do

              else if ((.NOT. flow_out(qk) .AND. &
                        bdry_face >= 101 .AND. bdry_face <= 200) .OR. &
                       bdry_face == 1) then

                call robin_bc(problem_dim, bdry_face, global_points_face(:, qk), gR)

                robin_term = -gR*phi1(ieqn, qk, i)

                ! RHS
                face_residual_p(ieqn, i) = face_residual_p(ieqn, i) &
                                           + integral_weighting(qk)*robin_term

              else

                write (io_err, "(A,I0)") "ERROR: Boundary not implemented, bdry_face ", bdry_face
                STOP - 1

              end if

            end do
          end do

! =======================================================================================

        else

! Interior face

          do qk = 1, no_quad_points
            interpolant_uh1(:, qk) = uh_face1(fe_basis_info, no_pdes, qk)
            interpolant_uh2(:, qk) = uh_face2(fe_basis_info, no_pdes, qk)

            do i = 1, no_pdes
              gradient_uh1(i, qk, 1:problem_dim) = grad_uh_face1(fe_basis_info, problem_dim, i, qk, 1)
              gradient_uh2(i, qk, 1:problem_dim) = grad_uh_face2(fe_basis_info, problem_dim, i, qk, 1)
            end do

            call compute_uh_glob_pt(adv_vel1(:, qk), &
                                    no_pdes_fluid, neighbours(1), &
                                    global_points_face(:, qk), problem_dim_fluid, mesh_data_fluid, soln_data_fluid)
            adv_vel1(problem_dim_fluid + 1, qk) = 0.0_db
            call compute_uh_glob_pt(adv_vel2(:, qk), &
                                    no_pdes_fluid, neighbours(2), &
                                    global_points_face(:, qk), problem_dim_fluid, mesh_data_fluid, soln_data_fluid)
            adv_vel2(problem_dim_fluid + 1, qk) = 0.0_db
            normal_vel(qk) = 0.5_db* &
                             dot_product(adv_vel1(1:problem_dim_fluid, qk) + adv_vel2(1:problem_dim_fluid, qk), &
                                         face_normals(:, qk))
            if (normal_vel(qk) > -tol) then
              flow_out(qk) = .TRUE.
            else
              flow_out(qk) = .FALSE.
            end if

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

! =======================================================================================
! Momentum Equations

          do ieqn = 1, no_pdes
            do qk = 1, no_quad_points

              if (flow_out(qk)) then
                upwind_val_int_face = interpolant_uh1(ieqn, qk)
              else
                upwind_val_int_face = interpolant_uh2(ieqn, qk)
              end if

! Loop over phi_i

              do i = 1, no_dofs_per_variable1(ieqn)

                adv_terms = -upwind_val_int_face*phi1(ieqn, qk, i)*normal_vel(qk)

                diff_terms = (1.0_db/Pe)*( &
                             -full_dispenal*(interpolant_uh1(ieqn, qk) - interpolant_uh2(ieqn, qk))*phi1(ieqn, qk, i) &
                             + 0.5_db* &
                         dot_product(gradient_uh1(ieqn, qk, :) + gradient_uh2(ieqn, qk, :), face_normals(:, qk))*phi1(ieqn, qk, i) &
                             + 0.5_db*dot_product(grad_phi1(ieqn, qk, :, i), face_normals(:, qk))* &
                             (interpolant_uh1(ieqn, qk) - interpolant_uh2(ieqn, qk)) &
                             )

                face_residual_p(ieqn, i) = face_residual_p(ieqn, i) + integral_weighting(qk)*( &
                                           adv_terms + diff_terms)

              end do

              do i = 1, no_dofs_per_variable2(ieqn)

                adv_terms = upwind_val_int_face*phi2(ieqn, qk, i)*normal_vel(qk)

                diff_terms = (1.0_db/Pe)*( &
                             full_dispenal*(interpolant_uh1(ieqn, qk) - interpolant_uh2(ieqn, qk))*phi2(ieqn, qk, i) &
                             - 0.5_db* &
                             dot_product( &
                             gradient_uh1(ieqn, qk, :) + gradient_uh2(ieqn, qk, :), face_normals(:, qk))*phi2(ieqn, qk, i) &
                             + 0.5_db*dot_product(grad_phi2(ieqn, qk, :, i), face_normals(:, qk))* &
                             (interpolant_uh1(ieqn, qk) - interpolant_uh2(ieqn, qk)) &
                             )

                face_residual_m(ieqn, i) = face_residual_m(ieqn, i) + integral_weighting(qk)*( &
                                           adv_terms + diff_terms)

              end do

            end do
          end do

! =======================================================================================

        end if

        deallocate (adv_vel1, adv_vel2)

      end associate

    end associate

  end subroutine adr_av_dg_face_nonlin_residual

!--------------------------------------------------------------------
! Jacobian for ADR
!--------------------------------------------------------------------
  subroutine adr_av_jac_element_mat(element_matrix, &
                                    mesh_data, soln_data, facet_data, fe_basis_info)
!--------------------------------------------------------------------
    use functions

    include 'assemble_jac_matrix_element.h'

! Local variables
    integer :: qk, i, j, ieqn, ivar
    real(db) :: adv_terms, diff_terms, reaction_terms, residual, diff_residual, &
                diff_residual_terms, ele_diam, coeff_AV, div_u, res_sgn
    real(db), dimension(facet_data%no_quad_points) :: eff_Dm
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points) :: interpolant_uh
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points, &
                        facet_data%problem_dim) :: gradient_uh
    real(db), dimension(facet_data%dim_soln_coeff, facet_data%no_quad_points, &
                        facet_data%problem_dim, maxval(facet_data%no_dofs_per_variable)) :: grad_phi
    real(db), dimension(facet_data%dim_soln_coeff, facet_data%no_quad_points, &
                        maxval(facet_data%no_dofs_per_variable)) :: phi
    real(db), dimension(:, :), allocatable :: adv_vel
    real(db), dimension(:, :, :), allocatable :: grad_adv_vel
    integer :: dim_soln_coeff_fluid, no_pdes_fluid, problem_dim_fluid

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

      associate ( &
        soln_data_fluid => fe_data_struct(1)%external_soln, &
        mesh_data_fluid => fe_data_struct(1)%external_mesh)
        dim_soln_coeff_fluid = get_dim_soln_coeff(soln_data_fluid)
        no_pdes_fluid = get_no_pdes(soln_data_fluid)
        problem_dim_fluid = get_problem_dim(mesh_data_fluid)
        allocate (adv_vel(no_pdes_fluid, facet_data%no_quad_points), &
                  grad_adv_vel(no_pdes_fluid, problem_dim_fluid, facet_data%no_quad_points))

        element_matrix = 0.0_db

        div_u = 0.0_db
        residual = 0.0_db
        ele_diam = compute_h_iso(mesh_data, element_number)
        coeff_AV = C_AV*(ele_diam**(2.0_db - B_AV))

        do qk = 1, no_quad_points
          interpolant_uh(:, qk) = uh_element(fe_basis_info, no_pdes, qk)

          do i = 1, no_pdes
            gradient_uh(i, qk, 1:problem_dim) = grad_uh_element(fe_basis_info, problem_dim, i, qk, 1)
          end do

          call compute_uh_gradient_uh_glob_pt(adv_vel(:, qk), grad_adv_vel(:, :, qk), &
                                              no_pdes_fluid, element_number, &
                                              global_points_ele(:, qk), problem_dim_fluid, mesh_data_fluid, soln_data_fluid)
          adv_vel(problem_dim_fluid + 1, qk) = 0.0_db

          call eval_tanh_smooth(problem_dim, element_number, global_points_ele(:, qk), &
                                Dm, eff_Dm(qk))

        end do

! Calculate Basis Functions

        do i = 1, dim_soln_coeff
          grad_phi(i, 1:no_quad_points, 1:problem_dim, 1:no_dofs_per_variable(i)) = fe_basis_info%basis_element &
                                        %deriv_basis_fns(i)%grad_data(1:no_quad_points, 1:problem_dim, 1:no_dofs_per_variable(i), 1)
          phi(i, 1:no_quad_points, 1:no_dofs_per_variable(i)) = fe_basis_info%basis_element%basis_fns(i) &
                                                                %fem_basis_fns(1:no_quad_points, 1:no_dofs_per_variable(i), 1)
        end do

        do qk = 1, no_quad_points

          do j = 1, problem_dim
            div_u = div_u + grad_adv_vel(j, j, qk)
          end do

! Loop over the equations

          do ieqn = 1, no_pdes

! Loop over phi_i

            do i = 1, no_dofs_per_variable(ieqn)

! Loop over the variables

              do ivar = 1, no_pdes

                residual = &
                  interpolant_uh(ivar, qk)*div_u &
                  + dot_product(adv_vel(1:problem_dim_fluid, qk), gradient_uh(ivar, qk, :)) &
                  + eff_Dm(qk)*interpolant_uh(ivar, qk)

                res_sgn = cal_sign(residual)
                residual = abs(residual)

! Loop over phi_j

                do j = 1, no_dofs_per_variable(ivar)

                  diff_residual = &
                    phi(ivar, qk, j)*div_u &
                    + dot_product(adv_vel(1:problem_dim_fluid, qk), grad_phi(ivar, qk, :, j)) &
                    + eff_Dm(qk)*phi(ivar, qk, j)

                  diff_residual_terms = &
                    coeff_AV*( &
                    res_sgn*diff_residual* &
                    dot_product(gradient_uh(ivar, qk, :), grad_phi(ieqn, qk, :, i)) &
                    + residual*dot_product( &
                    grad_phi(ivar, qk, :, j), grad_phi(ieqn, qk, :, i)) &
                    )

                  adv_terms = -phi(ivar, qk, j)* &
                              dot_product(adv_vel(1:problem_dim_fluid, qk), grad_phi(ieqn, qk, :, i))

                  diff_terms = (1.0_db/Pe)*( &
                               dot_product(grad_phi(ivar, qk, :, j), grad_phi(ieqn, qk, :, i)) &
                               )

                  reaction_terms = eff_Dm(qk)*phi(ivar, qk, j)*phi(ieqn, qk, i)

                  element_matrix(ieqn, ivar, i, j) = element_matrix(ieqn, ivar, i, j) &
                                                     + integral_weighting(qk)*(adv_terms + diff_terms + reaction_terms &
                                                                               + diff_residual_terms)

                end do

              end do

            end do

          end do

        end do

        deallocate (adv_vel, grad_adv_vel)

      end associate

    end associate

  end subroutine adr_av_jac_element_mat

!--------------------------------------------------------------------
! Face Jacobian for ADR
!--------------------------------------------------------------------
  subroutine adr_av_dg_jac_face_mat(face_matrix_pp, face_matrix_pm, &
                                    face_matrix_mp, face_matrix_mm, mesh_data, soln_data, &
                                    facet_data, fe_basis_info)
!--------------------------------------------------------------------
    use aptofem_fe_matrix_assembly

    include 'assemble_jac_matrix_int_bdry_face.h'

! Local variables

    integer :: qk, i, j, ieqn, ivar
    real(db) :: full_dispenal
    real(db), dimension(facet_data%no_pdes, facet_data%no_quad_points) :: uloc
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

    integer :: dim_soln_coeff_fluid, no_pdes_fluid, problem_dim_fluid
    logical, dimension(facet_data%no_quad_points) :: flow_out
    real(db) :: adv_terms, diff_terms, upwind_val_int_face, tol
    real(db), dimension(facet_data%no_quad_points) :: normal_vel
    real(db), dimension(:, :), allocatable :: &
      adv_vel1, adv_vel2

    tol = 1.0d-6

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

      associate ( &
        soln_data_fluid => fe_data_struct(1)%external_soln, &
        mesh_data_fluid => fe_data_struct(1)%external_mesh)
        dim_soln_coeff_fluid = get_dim_soln_coeff(soln_data_fluid)
        no_pdes_fluid = get_no_pdes(soln_data_fluid)
        problem_dim_fluid = get_problem_dim(mesh_data_fluid)
        allocate (adv_vel1(no_pdes_fluid, facet_data%no_quad_points), &
                  adv_vel2(no_pdes_fluid, facet_data%no_quad_points))

        face_matrix_pp = 0.0_db
        face_matrix_pm = 0.0_db
        face_matrix_mp = 0.0_db
        face_matrix_mm = 0.0_db

        dispenal_local = dispenal

        full_dispenal = interior_penalty_parameter*dispenal_local(1)*Pe

        if (bdry_face > 0) then

! ====================================================================
! Boundary Face

! Calculate value of analytical solution at quadrature points

          do qk = 1, no_quad_points
            interpolant_uh1(:, qk) = uh_face1(fe_basis_info, no_pdes, qk)
            do i = 1, no_pdes
              grad_uh1(i, qk, 1:problem_dim) = grad_uh_face1(fe_basis_info, problem_dim, i, qk, 1)
            end do

            do i = 1, no_pdes
              grad_phi1(i, 1:no_quad_points, 1:problem_dim, 1:no_dofs_per_variable1(i)) = fe_basis_info%basis_face1 &
                                                   %deriv_basis_fns(i)%grad_data(1:no_quad_points, :, 1:no_dofs_per_variable1(i), 1)
              phi1(i, 1:no_quad_points, 1:no_dofs_per_variable1(i)) = fe_basis_info%basis_face1%basis_fns(i) &
                                                                     %fem_basis_fns(1:no_quad_points, 1:no_dofs_per_variable1(i), 1)
            end do

            call compute_uh_glob_pt(adv_vel1(:, qk), &
                                    no_pdes_fluid, neighbours(1), &
                                    global_points_face(:, qk), problem_dim_fluid, mesh_data_fluid, soln_data_fluid)
            adv_vel1(problem_dim_fluid + 1, qk) = 0.0_db

            normal_vel(qk) = dot_product(adv_vel1(1:problem_dim_fluid, qk), face_normals(:, qk))
            if (normal_vel(qk) > -tol) then
              flow_out(qk) = .TRUE.
            else
              flow_out(qk) = .FALSE.
            end if

          end do

          do qk = 1, no_quad_points

! Loop over the equations

            do ieqn = 1, no_pdes

! Loop over phi_i

              do i = 1, no_dofs_per_variable1(ieqn)

! Loop over the variables

                do ivar = 1, no_pdes

! Loop over phi_j

                  do j = 1, no_dofs_per_variable1(ivar)

                    ! Neumann
                    if (flow_out(qk) .AND. &
                        bdry_face >= 101 .AND. bdry_face <= 200) then ! Neumann boundary

                      adv_terms = phi1(ivar, qk, j)*phi1(ieqn, qk, i)*normal_vel(qk)

                      diff_terms = 0.0_db

                      face_matrix_pp(ieqn, ivar, i, j) = face_matrix_pp(ieqn, ivar, i, j) &
                                                         + integral_weighting(qk)*(adv_terms + diff_terms)

                      ! Dirichlet
                    else if (bdry_face >= 2 .AND. bdry_face <= 100) then

                      adv_terms = 0.0_db

                      diff_terms = (1.0_db/Pe)*( &
                                   full_dispenal*phi1(ivar, qk, j)*phi1(ieqn, qk, i) &
                                   - dot_product(grad_phi1(ivar, qk, :, j), face_normals(:, qk))*phi1(ieqn, qk, i) &
                                   - dot_product(grad_phi1(ieqn, qk, :, i), face_normals(:, qk))*phi1(ivar, qk, j) &
                                   )

                      face_matrix_pp(ieqn, ivar, i, j) = face_matrix_pp(ieqn, ivar, i, j) &
                                                         + integral_weighting(qk)*(adv_terms + diff_terms)

                      ! Robin
                    else if ((.NOT. flow_out(qk) .AND. &
                              bdry_face >= 101 .AND. bdry_face <= 200) .OR. &
                             bdry_face == 1) then

                      adv_terms = 0.0_db

                      diff_terms = 0.0_db

                      face_matrix_pp(ieqn, ivar, i, j) = face_matrix_pp(ieqn, ivar, i, j) &
                                                         + integral_weighting(qk)*(adv_terms + diff_terms)

                    else

                      write (io_err, "(A,I0)") "ERROR: Boundary not implemented, bdry_face ", bdry_face
                      STOP - 1

                    end if

                  end do

                end do

              end do

            end do

          end do

        else

! ====================================================================
! Interior Face

          do qk = 1, no_quad_points
            interpolant_uh1(:, qk) = uh_face1(fe_basis_info, no_pdes, qk)
            interpolant_uh2(:, qk) = uh_face2(fe_basis_info, no_pdes, qk)
            do i = 1, no_pdes
              grad_uh1(i, qk, 1:problem_dim) = grad_uh_face1(fe_basis_info, problem_dim, i, qk, 1)
              grad_uh2(i, qk, 1:problem_dim) = grad_uh_face2(fe_basis_info, problem_dim, i, qk, 1)
            end do

            call compute_uh_glob_pt(adv_vel1(:, qk), &
                                    no_pdes_fluid, neighbours(1), &
                                    global_points_face(:, qk), problem_dim_fluid, mesh_data_fluid, soln_data_fluid)
            adv_vel1(problem_dim_fluid + 1, qk) = 0.0_db
            call compute_uh_glob_pt(adv_vel2(:, qk), &
                                    no_pdes_fluid, neighbours(2), &
                                    global_points_face(:, qk), problem_dim_fluid, mesh_data_fluid, soln_data_fluid)
            adv_vel2(problem_dim_fluid + 1, qk) = 0.0_db

            normal_vel(qk) = 0.5_db* &
                             dot_product(adv_vel1(1:problem_dim_fluid, qk) + adv_vel2(1:problem_dim_fluid, qk), &
                                         face_normals(:, qk))
            if (normal_vel(qk) > -tol) then
              flow_out(qk) = .TRUE.
            else
              flow_out(qk) = .FALSE.
            end if

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

! Loop over phi_i

                do i = 1, no_dofs_per_variable1(ieqn)

! u^+ v^+ w.r.t. ele1

                  do j = 1, no_dofs_per_variable1(ivar)

                    if (flow_out(qk)) then
                      upwind_val_int_face = phi1(ivar, qk, j)
                    else
                      upwind_val_int_face = 0.0_db
                    end if

                    adv_terms = upwind_val_int_face*phi1(ieqn, qk, i)*normal_vel(qk)

                    diff_terms = (1.0_db/Pe)*( &
                                 full_dispenal*phi1(ivar, qk, j)*phi1(ieqn, qk, i) &
                                 - 0.5_db* &
                                 dot_product(grad_phi1(ivar, qk, :, j), face_normals(:, qk))*phi1(ieqn, qk, i) &
                                 - 0.5_db*dot_product(grad_phi1(ieqn, qk, :, i), face_normals(:, qk))*phi1(ivar, qk, j) &
                                 )

                    face_matrix_pp(ieqn, ivar, i, j) = face_matrix_pp(ieqn, ivar, i, j) &
                                                       + integral_weighting(qk)*( &
                                                       adv_terms + diff_terms)

                  end do

! u^- v^+ w.r.t. ele1

! Loop over phi_j

                  do j = 1, no_dofs_per_variable2(ivar)

                    if (flow_out(qk)) then
                      upwind_val_int_face = 0.0_db
                    else
                      upwind_val_int_face = phi2(ivar, qk, j)
                    end if

                    adv_terms = upwind_val_int_face*phi1(ieqn, qk, i)*normal_vel(qk)

                    diff_terms = (1.0_db/Pe)*( &
                                 -full_dispenal*phi2(ivar, qk, j)*phi1(ieqn, qk, i) &
                                 - 0.5_db* &
                                 dot_product(grad_phi2(ivar, qk, :, j), face_normals(:, qk))*phi1(ieqn, qk, i) &
                                 + 0.5_db*dot_product(grad_phi1(ieqn, qk, :, i), face_normals(:, qk))*phi2(ivar, qk, j) &
                                 )

                    face_matrix_mp(ieqn, ivar, i, j) = face_matrix_mp(ieqn, ivar, i, j) &
                                                       + integral_weighting(qk)*( &
                                                       adv_terms + diff_terms)

                  end do

                end do

! Loop over phi_i

                do i = 1, no_dofs_per_variable2(ieqn)

! u^+ v^- w.r.t. ele1

! Loop over phi_j

                  do j = 1, no_dofs_per_variable1(ivar)

                    if (flow_out(qk)) then
                      upwind_val_int_face = phi1(ivar, qk, j)
                    else
                      upwind_val_int_face = 0.0_db
                    end if

                    adv_terms = -upwind_val_int_face*phi2(ieqn, qk, i)*normal_vel(qk)

                    diff_terms = (1.0_db/Pe)*( &
                                 -full_dispenal*phi1(ivar, qk, j)*phi2(ieqn, qk, i) &
                                 + 0.5_db* &
                                 dot_product(grad_phi1(ivar, qk, :, j), face_normals(:, qk))*phi2(ieqn, qk, i) &
                                 - 0.5_db*dot_product(grad_phi2(ieqn, qk, :, i), face_normals(:, qk))*phi1(ivar, qk, j) &
                                 )

                    face_matrix_pm(ieqn, ivar, i, j) = face_matrix_pm(ieqn, ivar, i, j) &
                                                       + integral_weighting(qk)*( &
                                                       adv_terms + diff_terms)

                  end do

! u^- v^- w.r.t. ele1

! Loop over phi_j

                  do j = 1, no_dofs_per_variable2(ivar)

                    if (flow_out(qk)) then
                      upwind_val_int_face = 0.0_db
                    else
                      upwind_val_int_face = phi2(ivar, qk, j)
                    end if

                    adv_terms = -upwind_val_int_face*phi2(ieqn, qk, i)*normal_vel(qk)

                    diff_terms = (1.0_db/Pe)*( &
                                 full_dispenal*phi2(ivar, qk, j)*phi2(ieqn, qk, i) &
                                 + 0.5_db* &
                                 dot_product(grad_phi2(ivar, qk, :, j), face_normals(:, qk))*phi2(ieqn, qk, i) &
                                 + 0.5_db*dot_product(grad_phi2(ieqn, qk, :, i), face_normals(:, qk))*phi2(ivar, qk, j) &
                                 )

                    face_matrix_mm(ieqn, ivar, i, j) = face_matrix_mm(ieqn, ivar, i, j) &
                                                       + integral_weighting(qk)*( &
                                                       adv_terms + diff_terms)

                  end do

                end do

              end do

            end do

          end do

! End of interior face
! ====================================================================
        end if

        deallocate (adv_vel1, adv_vel2)

      end associate

    end associate

  end subroutine adr_av_dg_jac_face_mat

end module adr_av_jac_matrix_residual
