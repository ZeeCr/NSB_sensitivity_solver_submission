! ------------------------------
! TODO: Read from CSV if flow profile ever changes
! ------------------------------
module input_data_storage

  use param

  implicit none

  integer, parameter :: no_data_pts_sa = 42

  real(db), dimension(no_data_pts_sa) :: time_values_sa, peak_vel_sa

contains

  subroutine initialise_input_data(t_scaling)

    implicit none

    real(db), intent(in), optional :: t_scaling

    time_values_sa = (/ &
                     0.0_db, &
                     0.09204_db, &
                     0.24187254_db, &
                     0.33523758_db, &
                     0.37437752_db, &
                     0.46723293_db, &
                     0.48863758_db, &
                     0.58098335_db, &
                     0.6019803_db, &
                     0.73061205_db, &
                     1.10998398_db, &
                     1.41647818_db, &
                     1.68607485_db, &
                     1.77413969_db, &
                     1.95291948_db, &
                     2.20427123_db, &
                     2.27398922_db, &
                     2.37989128_db, &
                     2.57589672_db, &
                     2.77302335_db, &
                     2.89737417_db, &
                     3.07533854_db, &
                     3.29009853_db, &
                     3.32373441_db, &
                     3.53808669_db, &
                     3.75254089_db, &
                     4.00307723_db, &
                     4.34483813_db, &
                     4.70494588_db, &
                     5.1740135_db, &
                     5.64379459_db, &
                     6.13243217_db, &
                     6.62056012_db, &
                     7.05456487_db, &
                     7.45136631_db, &
                     7.90137359_db, &
                     8.35209436_db, &
                     8.69528224_db, &
                     9.10900354_db, &
                     9.48643886_db, &
                     9.84675047_db, &
                     10.5_db &
                     /)

    peak_vel_sa = (/ &
                  0.97664803_db, &
                  1.170737778_db, &
                  1.348458238_db, &
                  1.450480642_db, &
                  1.549228058_db, &
                  1.634790676_db, &
                  1.74341962_db, &
                  1.812522448_db, &
                  1.907983562_db, &
                  1.983658994_db, &
                  2.0_db, &
                  1.993319928_db, &
                  1.960315508_db, &
                  1.891156116_db, &
                  1.838428232_db, &
                  1.798845552_db, &
                  1.719815944_db, &
                  1.644066982_db, &
                  1.56499778_db, &
                  1.522140112_db, &
                  1.459553324_db, &
                  1.380489778_db, &
                  1.324458622_db, &
                  1.245440328_db, &
                  1.176241342_db, &
                  1.110334314_db, &
                  1.044415974_db, &
                  1.011388928_db, &
                  0.9882321_db, &
                  0.988085036_db, &
                  1.010981676_db, &
                  1.06020832_db, &
                  1.092975178_db, &
                  1.125759002_db, &
                  1.122342606_db, &
                  1.089281624_db, &
                  1.079264344_db, &
                  1.092324704_db, &
                  1.052691118_db, &
                  1.00648493_db, &
                  0.989912016_db, &
                  0.97664803_db &
                  /)

    if (present(t_scaling)) then
      ! Scale time values by t_scaling
      time_values_sa = time_values_sa*t_scaling
    end if

  end subroutine initialise_input_data

end module input_data_storage
