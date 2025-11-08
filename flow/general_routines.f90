!--------------------------------------------------------------------
!! Note: Backflow stabilisation term is calculated in normal_vel_stab
!! semi-implicitly, regardless of bdf_semi_implicit flag
!--------------------------------------------------------------------
module general_routines

  implicit none

contains

  function cal_sym_stress(grad_u, grad_v, ieqn, problem_dim, no_pdes)

    use param
    implicit none

    integer, intent(in) :: problem_dim, no_pdes
    real(db) :: cal_sym_stress
    real(db), dimension(no_pdes, problem_dim), intent(in) :: grad_u, grad_v
    integer, intent(in) :: ieqn

    cal_sym_stress = 0.0_db

    if (ieqn < problem_dim + 1) then
      cal_sym_stress = 2.0_db*( &
                       dot_product(grad_u(ieqn, 1:problem_dim), grad_v(ieqn, 1:problem_dim)) &
                       + dot_product(grad_u(1:problem_dim, ieqn), grad_v(ieqn, 1:problem_dim)) &
                       )
    end if

  end function cal_sym_stress

  function cal_bdry_sym_stress_u(grad_u, normal, ieqn, problem_dim, no_pdes)

    use param
    implicit none

    integer, intent(in) :: problem_dim, no_pdes
    real(db) :: cal_bdry_sym_stress_u
    real(db), dimension(problem_dim), intent(in) :: normal
    real(db), dimension(no_pdes, problem_dim), intent(in) :: grad_u
    integer, intent(in) :: ieqn

    cal_bdry_sym_stress_u = 0.0_db

    if (ieqn < problem_dim + 1) then
      cal_bdry_sym_stress_u = dot_product( &
                              grad_u(ieqn, 1:problem_dim) + grad_u(1:problem_dim, ieqn), &
                              normal(:) &
                              )
    end if

  end function cal_bdry_sym_stress_u

  function cal_bdry_sym_stress_v(grad_v, u, u_D, normal, ieqn, problem_dim, no_pdes)

    use param
    implicit none

    integer, intent(in) :: problem_dim, no_pdes
    real(db) :: cal_bdry_sym_stress_v
    real(db), dimension(problem_dim), intent(in) :: normal, grad_v
    real(db), dimension(no_pdes), intent(in) :: u, u_D
    integer, intent(in) :: ieqn

    cal_bdry_sym_stress_v = 0.0_db

    if (ieqn < problem_dim + 1) then
      cal_bdry_sym_stress_v = dot_product( &
                              grad_v(1:problem_dim), normal(:))*(u(ieqn) - u_D(ieqn)) &
                              + dot_product(grad_v(1:problem_dim), u(1:problem_dim) - u_D(1:problem_dim))*normal(ieqn)
    end if

  end function cal_bdry_sym_stress_v

  function cal_sym_stress_deriv(grad_phi_u, grad_phi_v, &
                                ivar, ieqn, problem_dim, no_pdes)

    use param
    implicit none

    integer, intent(in) :: problem_dim, no_pdes
    real(db) :: cal_sym_stress_deriv
    real(db), dimension(no_pdes, problem_dim), intent(in) :: &
      grad_phi_u, grad_phi_v
    integer, intent(in) :: ivar, ieqn

    cal_sym_stress_deriv = 0.0_db

    if (ivar < problem_dim + 1 .AND. ieqn < problem_dim + 1) then
      cal_sym_stress_deriv = &
        2.0_db*grad_phi_u(ivar, ieqn)*grad_phi_v(ieqn, ivar)
    end if

    if (ieqn == ivar .AND. ivar < problem_dim + 1) then
      cal_sym_stress_deriv = cal_sym_stress_deriv &
                             + 2.0_db*dot_product(grad_phi_u(ivar, 1:problem_dim), grad_phi_v(ieqn, 1:problem_dim))
    end if

  end function cal_sym_stress_deriv

  function cal_stress_terms_bdry_face(grad_phi_u, grad_phi_v, phi_u, phi_v, &
                                      normal, penalisation, ivar, ieqn, problem_dim, no_pdes)

    use param
    implicit none

    integer, intent(in) :: problem_dim, no_pdes
    real(db) :: cal_stress_terms_bdry_face
    real(db), dimension(problem_dim), intent(in) :: &
      normal, grad_phi_u, grad_phi_v
    real(db), intent(in) :: phi_u, phi_v
    real(db), intent(in) :: penalisation
    integer, intent(in) :: ivar, ieqn

    cal_stress_terms_bdry_face = 0.0_db

    if (ieqn < problem_dim + 1 .AND. ivar < problem_dim + 1) then
      cal_stress_terms_bdry_face = &
        -grad_phi_u(ieqn)*normal(ivar)*phi_v &
        - grad_phi_v(ivar)*normal(ieqn)*phi_u
    end if

    if (ieqn == ivar .and. ivar < problem_dim + 1) then
      cal_stress_terms_bdry_face = cal_stress_terms_bdry_face &
                                   - ( &
                                   dot_product(grad_phi_u, normal)*phi_v &
                                   + dot_product(grad_phi_v, normal)*phi_u &
                                   ) &
                                   + penalisation*phi_u*phi_v
    end if

  end function cal_stress_terms_bdry_face

  function cal_stress_terms_int_face(grad_phi_u, grad_phi_v, phi_u, phi_v, &
                                     normal, penalisation, ivar, ieqn, problem_dim, no_pdes)

    use param
    implicit none

    integer, intent(in) :: problem_dim, no_pdes
    real(db) :: cal_stress_terms_int_face
    real(db), dimension(problem_dim), intent(in) :: &
      normal, grad_phi_u, grad_phi_v
    real(db), intent(in) :: phi_u, phi_v
    real(db), intent(in) :: penalisation
    integer, intent(in) :: ivar, ieqn

    cal_stress_terms_int_face = 0.0_db

    if (ieqn < problem_dim + 1 .AND. ivar < problem_dim + 1) then
      cal_stress_terms_int_face = -0.5_db*( &
                                  grad_phi_u(ieqn)*normal(ivar)*phi_v &
                                  + grad_phi_v(ivar)*normal(ieqn)*phi_u &
                                  )
    end if

    if (ieqn == ivar .and. ivar < problem_dim + 1) then
      cal_stress_terms_int_face = cal_stress_terms_int_face &
                                  - 0.5_db*( &
                                  dot_product(grad_phi_u, normal)*phi_v &
                                  + dot_product(grad_phi_v, normal)*phi_u &
                                  ) &
                                  + penalisation*phi_u*phi_v
    end if

  end function cal_stress_terms_int_face

  function cal_tang_stress_terms_bdry_face(grad_phi_u, grad_phi_v, phi_u, phi_v, &
                                           normal, surface_tensor, ivar, ieqn, problem_dim, no_pdes)

    use param
    implicit none

    integer, intent(in) :: problem_dim, no_pdes
    real(db) :: cal_tang_stress_terms_bdry_face
    real(db), dimension(problem_dim), intent(in) ::        normal
    real(db), dimension(problem_dim, problem_dim), intent(in) :: &
      surface_tensor
    real(db), dimension(problem_dim) :: grad_phi_u, grad_phi_v
    real(db), intent(in) :: phi_u, phi_v
    integer, intent(in) :: ivar, ieqn

! Local variables
    integer :: i, j

    cal_tang_stress_terms_bdry_face = 0.0_db

    if (ieqn < problem_dim + 1 .AND. ivar < problem_dim + 1) then
!n.(gradu+gradu^T).(I-nn).v
      cal_tang_stress_terms_bdry_face = -( &
                                        normal(ivar)*dot_product( &
                                        grad_phi_u(1:problem_dim), surface_tensor(1:problem_dim, ieqn) &
                                        )*phi_v &
                                        + dot_product(normal(1:problem_dim), grad_phi_u(1:problem_dim))* &
                                        surface_tensor(ivar, ieqn)*phi_v &
                                        )
!n.(gradv+gradv^T).(I-nn).u
      cal_tang_stress_terms_bdry_face = cal_tang_stress_terms_bdry_face - ( &
                                        normal(ieqn)*dot_product( &
                                        grad_phi_v(1:problem_dim), surface_tensor(1:problem_dim, ivar) &
                                        )*phi_u &
                                        + dot_product(normal(1:problem_dim), grad_phi_v(1:problem_dim))* &
                                        surface_tensor(ieqn, ivar)*phi_u &
                                        )
    end if

  end function cal_tang_stress_terms_bdry_face

  function cal_tangential_penalty_deriv(phi_u, phi_v, &
                                        surface_tensor, ivar, ieqn, problem_dim, no_pdes)

    use param
    implicit none

    integer, intent(in) :: problem_dim, no_pdes
    real(db) :: cal_tangential_penalty_deriv
    real(db), dimension(problem_dim, problem_dim), intent(in) :: &
      surface_tensor
    real(db), intent(in) :: phi_u, phi_v
    integer, intent(in) :: ivar, ieqn

    cal_tangential_penalty_deriv = 0.0_db

    if (ieqn < problem_dim + 1 .AND. ivar < problem_dim + 1) then
      cal_tangential_penalty_deriv = phi_u*phi_v* &
                                     dot_product( &
                                     surface_tensor(ivar, :), surface_tensor(ieqn, :) &
                                     )
    end if

  end function cal_tangential_penalty_deriv

  function cal_gradgradterm(grad_phi_u, grad_phi_v, ivar, ieqn, problem_dim, no_pdes)

    use param
    implicit none

    integer, intent(in) :: problem_dim, no_pdes
    real(db) :: cal_gradgradterm
    real(db), dimension(problem_dim), intent(in) :: grad_phi_u, grad_phi_v
    integer, intent(in) :: ivar, ieqn

    cal_gradgradterm = 0.0_db

    if (ieqn == ivar .and. ivar < problem_dim + 1) then
      cal_gradgradterm = dot_product(grad_phi_u, grad_phi_v)
    end if

  end function cal_gradgradterm

! --------------------------------------------------------------
! This routine calculates the components of the v . grad p term
! --------------------------------------------------------------
  function cal_gradterm(phi_v, grad_phi_p, ip, iv, problem_dim, no_pdes)

    use param
    implicit none

    integer, intent(in) :: problem_dim, no_pdes
    real(db) :: cal_gradterm
    real(db), dimension(problem_dim), intent(in) :: grad_phi_p
    real(db), intent(in) :: phi_v
    integer, intent(in) :: ip, iv

    cal_gradterm = 0.0_db

    if (ip == problem_dim + 1) then
      if (iv < problem_dim + 1) then
        cal_gradterm = phi_v*grad_phi_p(iv)
      end if
    end if

  end function cal_gradterm

  function cal_diff_terms_bdry_face(grad_phi_u, grad_phi_v, phi_u, phi_v, &
                                    normal, penalisation, ivar, ieqn, problem_dim, no_pdes)

    use param
    implicit none

    integer, intent(in) :: problem_dim, no_pdes
    real(db) :: cal_diff_terms_bdry_face
    real(db), dimension(problem_dim), intent(in) :: grad_phi_u, grad_phi_v, normal
    real(db), intent(in) :: phi_u, phi_v
    real(db), intent(in) :: penalisation
    integer, intent(in) :: ivar, ieqn

    cal_diff_terms_bdry_face = 0.0_db

    if (ieqn == ivar .and. ivar < problem_dim + 1) then
      cal_diff_terms_bdry_face = -dot_product(grad_phi_v, normal)*phi_u &
                                 - dot_product(grad_phi_u, normal)*phi_v + penalisation*phi_u*phi_v
    end if

  end function cal_diff_terms_bdry_face

  function cal_diff_terms_int_face(grad_phi_u, grad_phi_v, phi_u, phi_v, &
                                   normal, penalisation, ivar, ieqn, problem_dim, no_pdes)

    use param
    implicit none

    integer, intent(in) :: problem_dim, no_pdes
    real(db) :: cal_diff_terms_int_face
    real(db), dimension(problem_dim), intent(in) :: grad_phi_u, grad_phi_v, normal
    real(db), intent(in) :: phi_u, phi_v
    real(db), intent(in) :: penalisation
    integer, intent(in) :: ivar, ieqn

    cal_diff_terms_int_face = 0.0_db

    if (ieqn == ivar .and. ivar < problem_dim + 1) then
      cal_diff_terms_int_face = -0.5_db*(dot_product(grad_phi_v, normal)*phi_u &
                                         + dot_product(grad_phi_u, normal)*phi_v) + penalisation*phi_u*phi_v
    end if

  end function cal_diff_terms_int_face

  function cal_grad_terms_bdry_face(phi_p, phi_v, normal, ip, iv, problem_dim, no_pdes)

    use param
    implicit none

    integer, intent(in) :: problem_dim, no_pdes
    real(db) :: cal_grad_terms_bdry_face
    real(db), dimension(problem_dim), intent(in) :: normal
    real(db), intent(in) :: phi_p, phi_v
    integer, intent(in) :: ip, iv

    cal_grad_terms_bdry_face = 0.0_db

    if (ip == problem_dim + 1) then
      if (iv < problem_dim + 1) then
        cal_grad_terms_bdry_face = phi_v*phi_p*normal(iv)
      end if
    end if

  end function cal_grad_terms_bdry_face

! --------------------------------------------------------------
! This routine calculates the normal viscous stress,
! i.e. n_i ( grad_{j} u^i + grad^i u_{j} ) n_j,
! at specified number of points
! --------------------------------------------------------------
  subroutine compute_normal_viscous_stress(normal_stress, grad_u, normal, &
                                           problem_dim, no_pdes, no_points)

    use param
    implicit none

    integer, intent(in) :: problem_dim, no_pdes, no_points
    real(db), dimension(problem_dim, no_points), intent(in) :: normal
    real(db), dimension(no_pdes, no_points, problem_dim), intent(in) :: grad_u
    real(db), dimension(no_points), intent(out) :: normal_stress
    !< n_i(P^{ij})n_j evaluated at point k

! Local variables
    integer :: i, j, k

    normal_stress = 0.0_db

    do k = 1, no_points
      do j = 1, problem_dim
        do i = 1, problem_dim
          normal_stress(k) = normal_stress(k) + &
                             normal(i, k)*(grad_u(j, k, i) + grad_u(i, k, j))*normal(j, k)
        end do
      end do
    end do

  end subroutine compute_normal_viscous_stress

! --------------------------------------------------------------
! This routine calculates the viscous stress tensor,
! i.e. ( grad_{j} u^i + grad^i u_{j} ) e_i e^j,
! at specified number of points
! --------------------------------------------------------------
  subroutine compute_viscous_stress_tensor_u(stress_tensor, grad_u, &
                                             problem_dim, no_pdes, no_points)

    use param
    implicit none

    integer, intent(in) :: problem_dim, no_pdes, no_points
    real(db), dimension(no_pdes, no_points, problem_dim), intent(in) :: grad_u
    real(db), dimension(problem_dim, problem_dim, no_points), intent(out) :: &
      stress_tensor !< (P^{ij})e_i e_j evaluated at point k

! Local variables
    integer :: i, j, k

    stress_tensor = 0.0_db

    do k = 1, no_points
      do j = 1, problem_dim
        do i = 1, problem_dim
          stress_tensor(i, j, k) = grad_u(j, k, i) + grad_u(i, k, j)
        end do
      end do
    end do

  end subroutine compute_viscous_stress_tensor_u

! --------------------------------------------------------------
! This routine calculates the tangential stress component, i.e.
! n.P.(I-nn) at specified number of points
! --------------------------------------------------------------
  subroutine compute_stress_tangential_projection(tang_stress_vector, &
                                                  stress_tensor, surface_tensor, &
                                                  normal, problem_dim, no_points)

    use param
    implicit none

    integer, intent(in) :: problem_dim, no_points
    real(db), dimension(problem_dim, no_points), intent(in) :: normal
    real(db), dimension(problem_dim, problem_dim, no_points), intent(in) :: &
      stress_tensor, surface_tensor
    real(db), dimension(problem_dim, no_points), intent(out) :: &
      tang_stress_vector

! Local variables
    integer :: i, j, k
    real(db) :: component

    tang_stress_vector = 0.0_db

    do k = 1, no_points
      do i = 1, problem_dim
        component = 0.0_db
        do j = 1, problem_dim
          component = component &
                      + dot_product(normal(:, k), stress_tensor(:, j, k))*surface_tensor(j, i, k)
        end do
        tang_stress_vector(i, k) = component
      end do
    end do

  end subroutine compute_stress_tangential_projection

! --------------------------------------------------------------
! This routine calculates the tangential stress component, i.e.
! n.P.(I-nn) at specified number of points
! --------------------------------------------------------------
  subroutine compute_stress_tangential_projection_u(tang_stress_vector, &
                                                    stress_tensor, surface_tensor, normal, problem_dim, no_points)

    use param
    implicit none

    integer, intent(in) :: problem_dim, no_points
    real(db), dimension(problem_dim, no_points), intent(in) :: normal
    real(db), dimension(problem_dim, problem_dim, no_points), intent(in) :: &
      stress_tensor, surface_tensor
    real(db), dimension(problem_dim, no_points), intent(out) :: &
      tang_stress_vector

! Local variables
    integer :: qk

    tang_stress_vector = 0.0_db

    do qk = 1, no_points
      tang_stress_vector(:, qk) = matmul(matmul(normal(:, qk), stress_tensor(:, :, qk)), surface_tensor(:, :, qk))
    end do

  end subroutine compute_stress_tangential_projection_u

! --------------------------------------------------------------
! This routine calculates the tangential stress component, i.e.
! n.P.(I-nn) at specified number of points
! --------------------------------------------------------------
  subroutine compute_stress_tangential_projection_v(tang_stress_vector_phi, &
                                                    u, grad_phi, surface_tensor, normal, ieqn, problem_dim, no_pdes)

    use param
    implicit none

    integer, intent(in) :: ieqn, problem_dim, no_pdes
    real(db), dimension(problem_dim), intent(in) :: normal, grad_phi
    real(db), dimension(no_pdes), intent(in) :: u
    real(db), dimension(problem_dim, problem_dim), intent(in) :: &
      surface_tensor
    real(db), intent(out) :: tang_stress_vector_phi

    tang_stress_vector_phi = 0.0_db

    tang_stress_vector_phi = &
      normal(ieqn)*dot_product( &
      matmul(grad_phi(1:problem_dim), surface_tensor), u(1:problem_dim) &
      ) &
      + dot_product(grad_phi(1:problem_dim), normal)* &
      dot_product(surface_tensor(ieqn, 1:problem_dim), u(1:problem_dim))

  end subroutine compute_stress_tangential_projection_v

! --------------------------------------------------------------
! This routine calculates (u.(I-nn)).(v.(I-nn))
! --------------------------------------------------------------
  subroutine compute_tangential_penalty(tang_penalty, u, phi, &
                                        surface_tensor, ieqn, problem_dim, no_pdes)

    use param
    implicit none

    integer, intent(in) :: ieqn, problem_dim, no_pdes
    real(db), dimension(no_pdes), intent(in) :: u
    real(db), intent(in) :: phi
    real(db), dimension(problem_dim, problem_dim), intent(in) :: &
      surface_tensor
    real(db), intent(out) :: tang_penalty

    tang_penalty = 0.0_db

    tang_penalty = dot_product( &
                   matmul(u(1:problem_dim), surface_tensor), &
                   surface_tensor(ieqn, 1:problem_dim) &
                   )*phi

  end subroutine compute_tangential_penalty

! -------------------------------------------------------------
!> Backflow stabilisation: term is 0 on outflow Neumann
!! boundaries. Using BDF2, in the semi-implicit case, the
!! u.n term uses the extrapolated velocity
!!!TODO: implicit case
! -------------------------------------------------------------
  subroutine normal_vel_stab(bf_stab_term, &
                             vel, normal, &
                             problem_dim, no_pdes, bf_stab_param)

    use param
    implicit none

    integer, intent(in) :: problem_dim, no_pdes
    real(db) :: bf_stab_param
    real(db), dimension(problem_dim), intent(in) :: normal
    real(db), dimension(no_pdes), intent(in) :: &
      vel
    real(db), intent(out) :: &
      bf_stab_term

! Local variables
    real(db) :: velocity_normal, tol

    tol = 1.0d-9

    velocity_normal = dot_product( &
                      vel(1:problem_dim), normal(:))

! If velocity points into domain (u.n < 0) on Nuemann boundary
! then penalise, else this penalty term vanishes
    bf_stab_term = -bf_stab_param &
                   *(velocity_normal - abs(velocity_normal))/2.0_db

  end subroutine normal_vel_stab

! -------------------------------------------------------------
!> Calculate \vect{u} _n = (\vect{u} . \vect{n} ) \vect{n}
!! at given number of points
! -------------------------------------------------------------
  subroutine velocity_normal_vector(vel_normal_vector, vel, &
                                    normal, problem_dim, no_points)

    use param
    implicit none

    integer, intent(in) :: problem_dim, no_points
    real(db), dimension(problem_dim, no_points), intent(in) :: &
      vel, normal
    real(db), dimension(problem_dim, no_points), intent(out) :: &
      vel_normal_vector

! Local variables
    integer :: i, k
    real(db) :: vel_normal

    vel_normal = 0.0_db

    do k = 1, no_points
      vel_normal = dot_product(vel(:, k), normal(:, k))
      do i = 1, problem_dim
        vel_normal_vector(i, k) = vel_normal*normal(i, k)
      end do
    end do

  end subroutine velocity_normal_vector

! -------------------------------------------------------------
!> This routine defines the convective fluxes
!!
!! Author:
!!  Paul Houston
! -------------------------------------------------------------
  subroutine convective_fluxes(soln, fluxes, problem_dim, no_pdes)
! -------------------------------------------------------------

    use param
    implicit none

    integer, intent(in) :: problem_dim !< problem dimension
    integer, intent(in) :: no_pdes !< Number of variables in PDE system
    real(db), dimension(no_pdes), intent(in) :: soln
    real(db), intent(out), dimension(problem_dim, problem_dim) :: fluxes

! Local variables

    real(db), dimension(problem_dim) :: velocity
    integer :: i, j

    velocity = soln(1:problem_dim)

    do i = 1, problem_dim
      do j = 1, problem_dim
        fluxes(i, j) = velocity(i)*velocity(j)
      end do
    end do

  end subroutine convective_fluxes

! -------------------------------------------------------------
!> This routine defines the extrapolated convective fluxes
!> for semi-implicit BDF with Newton-Gregory backward interp.
! -------------------------------------------------------------
  subroutine extrap_adv_fluxes(soln_current, fluxes, &
                               problem_dim, no_pdes, soln_extrap)
! -------------------------------------------------------------

    use param
    implicit none

    integer, intent(in) :: problem_dim !< problem dimension
    integer, intent(in) :: no_pdes !< Number of variables in PDE system
    real(db), dimension(no_pdes), intent(in) :: soln_current
    real(db), dimension(no_pdes), intent(in) :: soln_extrap
    real(db), dimension(problem_dim, problem_dim), intent(out) :: fluxes

! Local variables
    real(db), dimension(problem_dim) :: velocity_current
    real(db), dimension(problem_dim) :: velocity_extrapolation
    integer :: i, j

    velocity_current(1:problem_dim) = soln_current(1:problem_dim)
    velocity_extrapolation(1:problem_dim) = soln_extrap(1:problem_dim)

    do i = 1, problem_dim
      do j = 1, problem_dim
        fluxes(i, j) = velocity_extrapolation(i)*velocity_current(j)
      end do
    end do

  end subroutine extrap_adv_fluxes

!  -------------------------------------------------------------
  subroutine jacobian_convective_fluxes(soln, fluxes_prime, problem_dim, no_pdes)
!  -------------------------------------------------------------
    use param

    implicit none

    integer, intent(in) :: problem_dim !< problem dimension
    integer, intent(in) :: no_pdes !< Number of variables in PDE system
    real(db), dimension(no_pdes), intent(in) :: soln
    real(db), dimension(problem_dim, problem_dim, problem_dim) :: &
      fluxes_prime
    real(db), dimension(problem_dim) :: velocity

    velocity = soln(1:problem_dim)

    if (problem_dim == 2) then

      fluxes_prime(1, 1, 1) = 2.0_db*velocity(1)
      fluxes_prime(1, 1, 2) = 0.0_db
      fluxes_prime(1, 2, 1) = velocity(2)
      fluxes_prime(1, 2, 2) = velocity(1)

      fluxes_prime(2, 1, 1) = velocity(2)
      fluxes_prime(2, 1, 2) = velocity(1)
      fluxes_prime(2, 2, 1) = 0.0_db
      fluxes_prime(2, 2, 2) = 2.0_db*velocity(2)

    else if (problem_dim == 3) then

      fluxes_prime(1, 1, 1) = 2.0_db*velocity(1)
      fluxes_prime(1, 1, 2) = 0.0_db
      fluxes_prime(1, 1, 3) = 0.0_db
      fluxes_prime(1, 2, 1) = velocity(2)
      fluxes_prime(1, 2, 2) = velocity(1)
      fluxes_prime(1, 2, 3) = 0.0_db
      fluxes_prime(1, 3, 1) = velocity(3)
      fluxes_prime(1, 3, 2) = 0.0_db
      fluxes_prime(1, 3, 3) = velocity(1)

      fluxes_prime(2, 1, 1) = velocity(2)
      fluxes_prime(2, 1, 2) = velocity(1)
      fluxes_prime(2, 1, 3) = 0.0_db
      fluxes_prime(2, 2, 1) = 0.0_db
      fluxes_prime(2, 2, 2) = 2.0_db*velocity(2)
      fluxes_prime(2, 2, 3) = 0.0_db
      fluxes_prime(2, 3, 1) = 0.0_db
      fluxes_prime(2, 3, 2) = velocity(3)
      fluxes_prime(2, 3, 3) = velocity(2)

      fluxes_prime(3, 1, 1) = velocity(3)
      fluxes_prime(3, 1, 2) = 0.0_db
      fluxes_prime(3, 1, 3) = velocity(1)
      fluxes_prime(3, 2, 1) = 0.0_db
      fluxes_prime(3, 2, 2) = velocity(3)
      fluxes_prime(3, 2, 3) = velocity(2)
      fluxes_prime(3, 3, 1) = 0.0_db
      fluxes_prime(3, 3, 2) = 0.0_db
      fluxes_prime(3, 3, 3) = 2.0_db*velocity(3)

    end if

  end subroutine jacobian_convective_fluxes

!  -------------------------------------------------------------
  subroutine jac_adv_flux_extrap(fluxes_prime, problem_dim, no_pdes, &
                                 soln)
!  -------------------------------------------------------------
    use param

    implicit none

    integer, intent(in) :: problem_dim !< problem dimension
    integer, intent(in) :: no_pdes !< Number of variables in PDE system
    real(db), dimension(no_pdes), intent(in) :: soln
    real(db), dimension(problem_dim, problem_dim, problem_dim), intent(out) :: &
      fluxes_prime

! Local variables
    real(db), dimension(problem_dim) :: velocity

    velocity(1:problem_dim) = soln(1:problem_dim)

    if (problem_dim == 2) then

      fluxes_prime(1, 1, 1) = velocity(1)
      fluxes_prime(1, 1, 2) = 0.0_db
      fluxes_prime(1, 2, 1) = 0.0_db
      fluxes_prime(1, 2, 2) = velocity(1)

      fluxes_prime(2, 1, 1) = velocity(2)
      fluxes_prime(2, 1, 2) = 0.0_db
      fluxes_prime(2, 2, 1) = 0.0_db
      fluxes_prime(2, 2, 2) = velocity(2)

    else if (problem_dim == 3) then

      fluxes_prime(1, 1, 1) = velocity(1)
      fluxes_prime(1, 1, 2) = 0.0_db
      fluxes_prime(1, 1, 3) = 0.0_db
      fluxes_prime(1, 2, 1) = 0.0_db
      fluxes_prime(1, 2, 2) = velocity(1)
      fluxes_prime(1, 2, 3) = 0.0_db
      fluxes_prime(1, 3, 1) = 0.0_db
      fluxes_prime(1, 3, 2) = 0.0_db
      fluxes_prime(1, 3, 3) = velocity(1)

      fluxes_prime(2, 1, 1) = velocity(2)
      fluxes_prime(2, 1, 2) = 0.0_db
      fluxes_prime(2, 1, 3) = 0.0_db
      fluxes_prime(2, 2, 1) = 0.0_db
      fluxes_prime(2, 2, 2) = velocity(2)
      fluxes_prime(2, 2, 3) = 0.0_db
      fluxes_prime(2, 3, 1) = 0.0_db
      fluxes_prime(2, 3, 2) = 0.0_db
      fluxes_prime(2, 3, 3) = velocity(2)

      fluxes_prime(3, 1, 1) = velocity(3)
      fluxes_prime(3, 1, 2) = 0.0_db
      fluxes_prime(3, 1, 3) = 0.0_db
      fluxes_prime(3, 2, 1) = 0.0_db
      fluxes_prime(3, 2, 2) = velocity(3)
      fluxes_prime(3, 2, 3) = 0.0_db
      fluxes_prime(3, 3, 1) = 0.0_db
      fluxes_prime(3, 3, 2) = 0.0_db
      fluxes_prime(3, 3, 3) = velocity(3)

    end if

  end subroutine jac_adv_flux_extrap

!  -------------------------------------------------------------
  subroutine jac_adv_flux_bdf(fluxes_prime, uh, problem_dim, no_pdes, &
                              bdf_order, bdf_semi_implicit, uh_previous_time_step)
!  -------------------------------------------------------------
    use param
    use bdf_timestepping

    implicit none

    logical, intent(in) :: bdf_semi_implicit
    integer, intent(in) :: problem_dim, bdf_order, no_pdes
    real(db), dimension(no_pdes), intent(in) :: uh
    real(db), dimension(no_pdes, bdf_order), intent(in) :: uh_previous_time_step
    real(db), dimension(problem_dim, problem_dim, problem_dim), intent(out) :: &
      fluxes_prime

! Local variables
    real(db), dimension(no_pdes) :: soln_extrap
    real(db), dimension(problem_dim) :: velocity

    velocity(1:problem_dim) = uh(1:problem_dim)

    if (bdf_semi_implicit) then

      call bdf_extrapolate(soln_extrap(:), &
                           uh_previous_time_step(:, :), no_pdes, bdf_order)

      call jac_adv_flux_extrap(fluxes_prime, problem_dim, no_pdes, &
                               soln_extrap)

    else

      call jacobian_convective_fluxes(uh, fluxes_prime, problem_dim, no_pdes)

    end if

  end subroutine jac_adv_flux_bdf

!  -------------------------------------------------------------
  subroutine compute_boundary_condition(boundary_condition, computed_soln, &
                                        analytical_soln, bdryno, problem_dim, no_pdes)
!  -------------------------------------------------------------
    use param

    implicit none

    integer, intent(in) :: problem_dim !< problem dimension
    integer, intent(in) :: no_pdes !< Number of variables in PDE system
    real(db), dimension(no_pdes), intent(out) :: boundary_condition
    real(db), dimension(no_pdes), intent(in) :: computed_soln, &
                                                analytical_soln
    integer, intent(in) :: bdryno

    if (bdryno <= 100) then

      boundary_condition = analytical_soln

    else if (bdryno <= 200) then

      boundary_condition = computed_soln

    else

      print *, 'compute_boundary_condition: Boundary number incorrect', bdryno
      stop

    end if

  end subroutine compute_boundary_condition
!  -------------------------------------------------------------
  subroutine compute_boundary_condition_bdf(boundary_condition, &
                                            uloc, computed_soln, face_normals, &
                                            problem_dim, no_pdes, bdry_face, global_points_face, current_time, &
                                            unloc, p_bc)
!  -------------------------------------------------------------
    use param
    use matrix_assembly_data_type

    implicit none

    integer, intent(in) :: problem_dim, no_pdes, bdry_face
    real(db), intent(in) :: current_time
    real(db), dimension(problem_dim), intent(in) :: &
      face_normals, global_points_face
    real(db), dimension(no_pdes), intent(out) :: boundary_condition
    real(db), dimension(no_pdes), intent(in) :: computed_soln
    real(db), dimension(no_pdes), intent(out) :: uloc
    real(db), intent(inout), optional :: p_bc
    real(db), dimension(problem_dim), intent(inout), optional :: unloc

    boundary_condition = 0.0_db

    if (bdry_face <= 100) then

      call anal_soln(uloc(:), global_points_face(:), &
                     problem_dim, no_pdes, bdry_face, current_time)

      boundary_condition = uloc

    else if (bdry_face <= 200) then

      if (present(unloc)) then
        call neumann_bc(unloc(:), p_bc, &
                        global_points_face(:), problem_dim, bdry_face)
      end if

      boundary_condition = computed_soln

    else if (bdry_face > 300 .AND. bdry_face <= 500) then

      boundary_condition(1:problem_dim) = &
        dot_product(computed_soln(1:problem_dim), face_normals) &
        *face_normals(1:problem_dim)

    else

      write (io_err, *) "compute_boundary_condition_bdf error"
      write (io_err, *) "Incorrect boundary number"
      stop

    end if

  end subroutine compute_boundary_condition_bdf

!  -------------------------------------------------------------
  subroutine compute_jac_boundary_condition(boundary_jacobian, &
                                            computed_soln, bdryno, normal, problem_dim, no_pdes)
!  -------------------------------------------------------------
    use param

    implicit none

    integer, intent(in) :: problem_dim !< problem dimension
    integer, intent(in) :: no_pdes !< Number of variables in PDE system
    real(db), dimension(problem_dim, problem_dim), intent(out) :: &
      boundary_jacobian
    real(db), dimension(no_pdes), intent(in) :: computed_soln
    real(db), dimension(problem_dim), intent(in) :: normal
    integer, intent(in) :: bdryno
    integer :: iv, j

    if (bdryno <= 100) then

      boundary_jacobian = 0.0_db

    else if (bdryno <= 200) then

      boundary_jacobian = 0.0_db
      do iv = 1, problem_dim
        boundary_jacobian(iv, iv) = 1.0_db
      end do

    else if (bdryno >= 301 .AND. bdryno <= 500) then

      boundary_jacobian = 0.0_db
      do iv = 1, problem_dim
        do j = 1, problem_dim
          boundary_jacobian(iv, j) = &
            normal(iv)*normal(j)
        end do
      end do

    else

      print *, 'compute_jacobian_boundary_condition:'
      print *, 'Boundary number incorrect', bdryno
      stop

    end if

  end subroutine compute_jac_boundary_condition

! --------------------------------------------------------------
  subroutine lax_friedrichs(nflxsoln, u1, u2, normal, problem_dim, no_pdes)
! --------------------------------------------------------------
    use param

    implicit none

    integer, intent(in) :: problem_dim !< problem dimension
    integer, intent(in) :: no_pdes !< Number of variables in PDE system
    real(db), dimension(problem_dim), intent(out) :: nflxsoln
    real(db), dimension(no_pdes), intent(in) :: u1, u2
    real(db), dimension(problem_dim), intent(in) :: normal
    real(db), dimension(problem_dim, problem_dim) :: fluxes1, fluxes2
    real(db) :: alpha
    integer :: i

    call convective_fluxes(u1, fluxes1, problem_dim, no_pdes)
    call convective_fluxes(u2, fluxes2, problem_dim, no_pdes)

    alpha = cal_alpha(u1, u2, normal, problem_dim, no_pdes)

    nflxsoln = 0.0_db

    do i = 1, problem_dim
      nflxsoln = nflxsoln + (fluxes1(i, :) + fluxes2(i, :))*normal(i)
    end do

    nflxsoln = 0.5_db*(nflxsoln - alpha*(u2(1:problem_dim) - u1(1:problem_dim)))

  end subroutine lax_friedrichs

! --------------------------------------------------------------
  subroutine lax_friedrichs_bdf_bdry(nflxsoln, u1, u2, normal, &
                                     problem_dim, no_pdes, bdry_face, &
                                     bdf_order, bdf_semi_implicit, uh_previous_time_step1)
! --------------------------------------------------------------
!< This routine calculates the Lax Friedrichs flux for BDF
!< time discretisation on a boundary face
! --------------------------------------------------------------
    use param
    use matrix_assembly_data_type
    use bdf_timestepping

    implicit none

    logical, intent(in) :: bdf_semi_implicit
    integer, intent(in) :: problem_dim, no_pdes, bdry_face, bdf_order
    real(db), dimension(problem_dim), intent(out) :: nflxsoln
    real(db), dimension(no_pdes), intent(in) :: u1, u2
    real(db), dimension(problem_dim), intent(in) :: normal
    real(db), dimension(no_pdes, bdf_order) :: uh_previous_time_step1

! Local variables
    real(db), dimension(no_pdes) :: u1_extrap
    real(db), dimension(problem_dim, problem_dim) :: fluxes1, fluxes2
    real(db) :: alpha
    integer :: i

    nflxsoln = 0.0_db

    if (bdf_semi_implicit) then

      call bdf_extrapolate(u1_extrap(:), &
                           uh_previous_time_step1(:, :), no_pdes, bdf_order)

      if (bdry_face > 0 .AND. bdry_face <= 100) then

        ! No need to extrapolate u_{\Gamma} on Dirichlet boundaries
        alpha = cal_alpha(u1_extrap, u2, &
                          normal, problem_dim, no_pdes)
        call extrap_adv_fluxes(u1_extrap, fluxes1, &
                               problem_dim, no_pdes, u1)
        call extrap_adv_fluxes(u2, fluxes2, &
                               problem_dim, no_pdes, u2)

      else if (bdry_face > 100 .AND. bdry_face <= 200) then

        alpha = cal_alpha(u1_extrap, u1_extrap, &
                          normal, problem_dim, no_pdes)
        call extrap_adv_fluxes(u1_extrap, fluxes1, &
                               problem_dim, no_pdes, u1)
        fluxes2 = fluxes1

      else if (bdry_face > 300 .AND. bdry_face <= 500) then

        alpha = cal_alpha(u1_extrap, u1_extrap, &
                          normal, problem_dim, no_pdes)
        call extrap_adv_fluxes(u1_extrap, fluxes1, &
                               problem_dim, no_pdes, u1)
        call extrap_adv_fluxes(u1_extrap, fluxes2, &
                               problem_dim, no_pdes, u2)

      end if

      do i = 1, problem_dim
        nflxsoln = nflxsoln + (fluxes1(:, i) + fluxes2(:, i))*normal(i)
      end do

      nflxsoln = 0.5_db*(nflxsoln + alpha*(u1(1:problem_dim) - u2(1:problem_dim)))

    else

      call convective_fluxes(u1, fluxes1, problem_dim, no_pdes)
      call convective_fluxes(u2, fluxes2, problem_dim, no_pdes)

      alpha = cal_alpha(u1, u2, normal, problem_dim, no_pdes)

      do i = 1, problem_dim
        nflxsoln = nflxsoln + (fluxes1(i, :) + fluxes2(i, :))*normal(i)
      end do

      nflxsoln = 0.5_db*(nflxsoln - alpha*(u2(1:problem_dim) - u1(1:problem_dim)))

    end if

  end subroutine lax_friedrichs_bdf_bdry

! --------------------------------------------------------------
  subroutine lax_friedrichs_bdf_int(nflxsoln, u1, u2, normal, &
                                    problem_dim, no_pdes, bdry_face, &
                                    bdf_order, bdf_semi_implicit, uh_previous_time_step1, &
                                    uh_previous_time_step2)
! --------------------------------------------------------------
!< This routine calculates the Lax Friedrichs flux for BDF
!< time discretisation on an interior face
! --------------------------------------------------------------
    use param
    use matrix_assembly_data_type
    use bdf_timestepping

    implicit none

    logical, intent(in) :: bdf_semi_implicit
    integer, intent(in) :: problem_dim, no_pdes, bdry_face, bdf_order
    real(db), dimension(problem_dim), intent(out) :: nflxsoln
    real(db), dimension(no_pdes), intent(in) :: u1, u2
    real(db), dimension(problem_dim), intent(in) :: normal
    real(db), dimension(no_pdes, bdf_order) :: &
      uh_previous_time_step1, uh_previous_time_step2

! Local variables
    real(db), dimension(problem_dim) :: u1_normal_vector
    real(db), dimension(no_pdes) :: soln_extrap1, soln_extrap2
    real(db), dimension(problem_dim, problem_dim) :: fluxes1, fluxes2
    real(db) :: alpha
    integer :: i

    nflxsoln = 0.0_db

    if (bdf_semi_implicit) then

      call bdf_extrapolate(soln_extrap1(:), &
                           uh_previous_time_step1(:, :), no_pdes, bdf_order)
      call bdf_extrapolate(soln_extrap2(:), &
                           uh_previous_time_step2(:, :), no_pdes, bdf_order)

      alpha = cal_alpha(soln_extrap1, soln_extrap2, &
                        normal, problem_dim, no_pdes)
      call extrap_adv_fluxes(soln_extrap1, fluxes1, &
                             problem_dim, no_pdes, u1)
      call extrap_adv_fluxes(soln_extrap2, fluxes2, &
                             problem_dim, no_pdes, u2)

      do i = 1, problem_dim
        nflxsoln = nflxsoln + (fluxes1(:, i) + fluxes2(:, i))*normal(i)
      end do

      nflxsoln = 0.5_db*(nflxsoln + alpha*(u1(1:problem_dim) - u2(1:problem_dim)))

    else

      call convective_fluxes(u1, fluxes1, problem_dim, no_pdes)
      call convective_fluxes(u2, fluxes2, problem_dim, no_pdes)

      alpha = cal_alpha(u1, u2, normal, problem_dim, no_pdes)

      do i = 1, problem_dim
        nflxsoln = nflxsoln + (fluxes1(i, :) + fluxes2(i, :))*normal(i)
      end do

      nflxsoln = 0.5_db*(nflxsoln - alpha*(u2(1:problem_dim) - u1(1:problem_dim)))

    end if

  end subroutine lax_friedrichs_bdf_int

! --------------------------------------------------------------
  subroutine lax_friedrichs_jac_bdf_bdry( &
    boundary_jacobian, fluxes_prime1, fluxes_prime2, alpha, &
    u1, u2, face_normals, problem_dim, no_pdes, bdry_face, &
    bdf_order, bdf_semi_implicit, uh_previous_time_step1)
! --------------------------------------------------------------
!< This routine calculates the Jacobian of the Lax Friedrichs
!< flux for BDF time discretisation on a boundary face
! --------------------------------------------------------------
    use param
    use matrix_assembly_data_type
    use bdf_timestepping

    implicit none

    logical, intent(in) :: bdf_semi_implicit
    integer, intent(in) :: problem_dim, no_pdes, bdry_face, bdf_order
    real(db), intent(out) :: alpha
    real(db), dimension(no_pdes), intent(in) :: u1, u2
    real(db), dimension(problem_dim), intent(in) :: face_normals
    real(db), dimension(no_pdes, bdf_order) :: uh_previous_time_step1
    real(db), dimension(problem_dim, problem_dim), intent(out) :: boundary_jacobian
    real(db), dimension(problem_dim, problem_dim, problem_dim), intent(out) :: &
      fluxes_prime1, fluxes_prime2

! Local variables
    real(db), dimension(no_pdes) :: u1_extrap, u2_extrap

    u1_extrap = 0.0_db
    u2_extrap = 0.0_db

    if (bdf_semi_implicit) then

      call bdf_extrapolate(u1_extrap, &
                           uh_previous_time_step1, no_pdes, bdf_order)

      if (bdry_face > 0 .AND. bdry_face <= 100) then

        ! Do not need to extrapolate u_{\Gamma} on Dirichlet boundaries
        alpha = cal_alpha(u1_extrap, u2, &
                          face_normals, problem_dim, no_pdes)
        call jac_adv_flux_extrap(fluxes_prime1, problem_dim, no_pdes, &
                                 u1_extrap)
        call jac_adv_flux_extrap(fluxes_prime2, problem_dim, no_pdes, &
                                 u2)

      else if (bdry_face > 100 .AND. bdry_face <= 200) then

        alpha = cal_alpha(u1_extrap, u1_extrap, &
                          face_normals, problem_dim, no_pdes)
        call jac_adv_flux_extrap(fluxes_prime1, problem_dim, no_pdes, &
                                 u1_extrap)
        fluxes_prime2 = fluxes_prime1

      else if (bdry_face > 300 .AND. bdry_face <= 500) then

        alpha = cal_alpha(u1_extrap, u1_extrap, &
                          face_normals, problem_dim, no_pdes)
        call jac_adv_flux_extrap(fluxes_prime1, problem_dim, no_pdes, &
                                 u1_extrap)
        u2_extrap(1:problem_dim) = &
          dot_product(u1_extrap(1:problem_dim), face_normals)*face_normals
        call jac_adv_flux_extrap(fluxes_prime2, problem_dim, no_pdes, &
                                 u2_extrap)

      else

        write (io_err, *) "lax_friedrichs_jac_bdf_bdry error"
        write (io_err, *) "Incorrect boundary number"
        stop

      end if

      call jacobian_lf_flux_bdry(boundary_jacobian, &
                                 bdry_face, fluxes_prime2, alpha, face_normals, &
                                 u1_extrap, problem_dim, no_pdes)

    else

      alpha = cal_alpha(u1, u2, &
                        face_normals, problem_dim, no_pdes)
      call jacobian_convective_fluxes(u1, &
                                      fluxes_prime1, problem_dim, no_pdes)
      call jacobian_convective_fluxes(u2, &
                                      fluxes_prime2, problem_dim, no_pdes)
      call jacobian_lf_flux_bdry(boundary_jacobian, &
                                 bdry_face, fluxes_prime2, alpha, face_normals, &
                                 u1, problem_dim, no_pdes)

    end if

  end subroutine lax_friedrichs_jac_bdf_bdry

  function cal_alpha(u1, u2, normal, problem_dim, no_pdes)

    use param

    implicit none

    real(db) :: cal_alpha
    integer, intent(in) :: problem_dim, no_pdes
    real(db), dimension(problem_dim), intent(in) :: normal
    real(db), dimension(no_pdes), intent(in) :: u1, u2

    cal_alpha = 2.0_db*max(abs(dot_product(u1(1:problem_dim), normal)), &
                           abs(dot_product(u2(1:problem_dim), normal)))

  end function cal_alpha
!  -------------------------------------------------------------
  subroutine jacobian_lf_flux_bdry(bc_jac, bdryno, fluxes_prime, alpha, normal, &
                                   soln, problem_dim, no_pdes)
!  -------------------------------------------------------------
    use param

    implicit none

    integer, intent(in) :: problem_dim !< problem dimension
    integer, intent(in) :: no_pdes
    real(db), dimension(problem_dim, problem_dim), intent(out) :: bc_jac
    real(db), dimension(problem_dim, problem_dim, problem_dim), intent(in) :: fluxes_prime
    real(db), dimension(problem_dim), intent(in) :: normal
    real(db), dimension(no_pdes), intent(in) :: soln
    real(db), intent(in) :: alpha
    integer, intent(in) :: bdryno
    real(db), dimension(problem_dim, problem_dim) :: var_jac

    real(db), dimension(problem_dim, problem_dim) :: temp_mat
    integer :: i

    call compute_jac_boundary_condition(var_jac, &
                                        soln, bdryno, normal, problem_dim, no_pdes)

    bc_jac = 0.0_db

    temp_mat = 0.0_db
    do i = 1, problem_dim
      temp_mat = temp_mat + fluxes_prime(i, :, :)*normal(i)
    end do

    bc_jac = 0.5_db*(matmul(temp_mat, var_jac) - alpha*var_jac)

  end subroutine jacobian_lf_flux_bdry

end module general_routines
