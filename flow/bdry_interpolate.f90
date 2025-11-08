subroutine poiseuille_interpolate(scheme_coupled_data)

  use param
  use linear_algebra
  use coupled_data_type
  use input_data_storage

  implicit none

  real(db) :: tol = 1.0d-9

  real(db) :: peak_V

  real(db) :: shifted_time, interp_point, interp_value
  integer :: bdry_no, bdry_counter, no_data_pts

  type(coupled_user_data), intent(inout) :: scheme_coupled_data

  peak_V = 0.0_db

  associate (t => scheme_coupled_data%current_time)

    do bdry_counter = 1, size(scheme_coupled_data%bdry_map)

      bdry_no = scheme_coupled_data%bdry_map(bdry_counter)

      if (bdry_no > 1 .AND. bdry_no < 100) then ! inlet

        call calc_linear_interp(no_data_pts_sa, t, time_values_sa, peak_vel_sa, interp_value)

        ! Interpolated peak_V here
        peak_V = interp_value

      else

        write (io_err, *) "Error: bdry_interpolate.f90"
        write (io_err, *) "Unexpected bdry_no in bdry_map"
        stop -2

      end if

      scheme_coupled_data%peak_v(bdry_counter) = peak_V

    end do

  end associate

end subroutine poiseuille_interpolate

! Loop through array ,y_arr, to find closest element smaller than value
subroutine calc_linear_interp(size_of_arr, x_val, x_arr, y_arr, interp_val)

  use param

  implicit none

  integer, intent(in) :: size_of_arr
  real(db), intent(in) :: x_val
  real(db), dimension(size_of_arr), intent(in) :: x_arr, y_arr

  real(db), intent(out) :: interp_val

  ! Local variables
  logical :: found_index
  integer :: loop_index, arr_index
  real(db) :: tol, shifted_val, period_val, delta_x, linear_grad

  tol = 1.0d-9

  ! Debug - set to TRUE if interp finds an index to interpolate against
  found_index = .FALSE.

  ! Periodic array, final entry == time per periodic
  period_val = x_arr(size_of_arr)

  shifted_val = x_val - floor(x_val/period_val)*period_val

  do loop_index = 1, size_of_arr - 1

    if (shifted_val > x_arr(loop_index) - tol .AND. shifted_val < x_arr(loop_index + 1) + tol) then

      found_index = .TRUE.
      arr_index = loop_index
      exit

    end if

  end do

  if (.NOT. found_index) then
    print *, "ERROR: calc_linear_interp"
    print *, "Unable to find array index for point ", x_val
    call exit(-1)
  else

    delta_x = shifted_val - x_arr(arr_index)
    linear_grad = (y_arr(arr_index + 1) - y_arr(arr_index))/(x_arr(arr_index + 1) - x_arr(arr_index))

    interp_val = y_arr(arr_index) + delta_x*linear_grad

  end if

end subroutine calc_linear_interp
