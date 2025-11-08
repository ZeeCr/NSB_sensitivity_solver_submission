module boundary_data_storage

  use param

  type :: boundary_storage
    integer :: bdry_no
    real(db) :: radius
    real(db), dimension(3) :: centre, normal
  end type

  type :: cavity_storage
    integer :: cavity_no
    real(db) :: minor_axis, major_axis
    real(db), dimension(3) :: centre, orientation_normal
    real(db), dimension(3, 3) :: basis_map ! Maps to ellipsoid in standard Cartesian basis at (0,0,0)
    ! with major axis oriented along z-axis
  end type

  type :: outlet_storage
    integer :: outlet_no
    real(db) :: radius, fillet_radius, total_radius, vein_length
    real(db), dimension(3) :: centre, interior_face_centre, outward_unit_normal
  end type

  integer :: no_bdries, no_cavities, no_outlets
  integer, dimension(:), allocatable :: element_cavity_no
  integer, dimension(:), allocatable :: element_outlet_no
  real(db), parameter :: smoothing_freq = 3.68055 ! y0 = 0.90, alpha = 0.4
  type(boundary_storage), dimension(:), allocatable :: bdry_storage
  type(cavity_storage), dimension(:), allocatable :: cvty_storage
  type(outlet_storage), dimension(:), allocatable :: otlt_storage

contains

  subroutine setup_bdry_storage()

    implicit none

    character(len=256) line
    integer :: io_code, line_no, ipos1, ipos2, counter

    no_bdries = 0

    open (unit=90, file='../data/mesh_info.csv', iostat=io_code, status='old')

    if (io_code /= 0) then
      write (io_err, *) 'Reading file mesh_info.csv failed: ', io_code
      stop
    end if

    read (90, *)

    do
      read (90, '(A256)', iostat=io_code) line
      if (io_code /= 0) then
        exit
      else
        no_bdries = no_bdries + 1
      end if
    end do

    allocate (bdry_storage(no_bdries))

    rewind (unit=90)

    read (90, *)

    do line_no = 1, no_bdries
      read (90, '(A256)') line

      ipos1 = 1
      ipos2 = 2
      counter = 0
      do while (ipos2 > 1)
        ipos2 = index(line(ipos1:), ',')

        if (counter == 0) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) bdry_storage(line_no)%bdry_no
        else if (counter == 1) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) bdry_storage(line_no)%radius
        else if (counter == 2) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) bdry_storage(line_no)%centre(1)
        else if (counter == 3) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) bdry_storage(line_no)%centre(2)
        else if (counter == 4) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) bdry_storage(line_no)%centre(3)
        else if (counter == 5) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) bdry_storage(line_no)%normal(1)
        else if (counter == 6) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) bdry_storage(line_no)%normal(2)
        else if (ipos2 == 0) then
          read (line(ipos1:), *) bdry_storage(line_no)%normal(3)
        else
          write (io_err, *) "Error: setup_bdry_storage"
          stop -10
        end if

        ipos1 = ipos1 + ipos2
        counter = counter + 1

      end do

    end do

    close (unit=90)

  end subroutine setup_bdry_storage

  subroutine delete_bdry_storage()

    implicit none

    deallocate (bdry_storage)

  end subroutine delete_bdry_storage

  subroutine setup_cvty_storage()

    implicit none

    character(len=256) line
    integer :: io_code, line_no, ipos1, ipos2, counter

    no_cavities = 0

    open (unit=91, file='../data/cavity_info.csv', iostat=io_code, status='old')

    if (io_code /= 0) then
      write (io_err, *) 'Reading file cavity_info.csv failed: ', io_code
      STOP - 1
    end if

    read (91, *)

    do
      read (91, '(A256)', iostat=io_code) line
      if (io_code /= 0) then
        exit
      else
        no_cavities = no_cavities + 1
      end if
    end do

    allocate (cvty_storage(no_cavities))

    rewind (unit=91)

    read (91, *)

    do line_no = 1, no_cavities
      read (91, '(A256)') line

      ipos1 = 1
      ipos2 = 2
      counter = 0
      do while (ipos2 > 1)
        ipos2 = index(line(ipos1:), ',')

        if (counter == 0) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) cvty_storage(line_no)%cavity_no
        else if (counter == 1) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) cvty_storage(line_no)%minor_axis
        else if (counter == 2) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) cvty_storage(line_no)%major_axis
        else if (counter == 3) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) cvty_storage(line_no)%centre(1)
        else if (counter == 4) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) cvty_storage(line_no)%centre(2)
        else if (counter == 5) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) cvty_storage(line_no)%centre(3)
        else if (counter == 6) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) cvty_storage(line_no)%orientation_normal(1)
        else if (counter == 7) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) cvty_storage(line_no)%orientation_normal(2)
        else if (ipos2 == 0) then
          read (line(ipos1:), *) cvty_storage(line_no)%orientation_normal(3)
        else
          write (io_err, *) "Error: setup_cvty_storage"
          stop - 10
        end if

        ipos1 = ipos1 + ipos2
        counter = counter + 1

      end do

    end do

    close (unit=91)

    call setup_basis_map()

  end subroutine setup_cvty_storage

  subroutine delete_cvty_storage()

    implicit none

    deallocate (cvty_storage)

    if (allocated(element_cavity_no)) then
      element_cavity_no = 0
      deallocate (element_cavity_no)
    end if

  end subroutine delete_cvty_storage

  subroutine setup_otlt_storage()

    implicit none

    character(len=256) line
    integer :: io_code, line_no, ipos1, ipos2, counter

    no_outlets = 0

    open (unit=92, file='../data/outlet_info.csv', iostat=io_code, status='old')

    if (io_code /= 0) then
      write (io_err, *) 'Reading file outlet_info.csv failed: ', io_code
      write (io_err, *) 'Assuming no outlet file exists'
      close (unit=92)
      return
    end if

    read (92, *)

    do
      read (92, '(A256)', iostat=io_code) line
      if (io_code /= 0) then
        exit
      else
        no_outlets = no_outlets + 1
      end if
    end do

    allocate (otlt_storage(no_outlets))

    rewind (unit=92)

    read (92, *)

    do line_no = 1, no_outlets
      read (92, '(A256)') line

      ipos1 = 1
      ipos2 = 2
      counter = 0
      do while (ipos2 > 1)
        ipos2 = index(line(ipos1:), ',')

        if (counter == 0) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) otlt_storage(line_no)%outlet_no
        else if (counter == 1) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) otlt_storage(line_no)%radius
        else if (counter == 2) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) otlt_storage(line_no)%centre(1)
        else if (counter == 3) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) otlt_storage(line_no)%centre(2)
        else if (counter == 4) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) otlt_storage(line_no)%centre(3)
        else if (counter == 5) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) otlt_storage(line_no)%outward_unit_normal(1)
        else if (counter == 6) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) otlt_storage(line_no)%outward_unit_normal(2)
        else if (counter == 7) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) otlt_storage(line_no)%outward_unit_normal(3)
        else if (counter == 8) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) otlt_storage(line_no)%vein_length
        else if (ipos2 == 0) then
          read (line(ipos1:), *) otlt_storage(line_no)%fillet_radius
        else
          write (io_err, *) "Error: setup_outlet_storage"
          stop - 10
        end if

        ipos1 = ipos1 + ipos2
        counter = counter + 1

      end do

      otlt_storage(line_no)%total_radius = otlt_storage(line_no)%radius &
                                           + otlt_storage(line_no)%fillet_radius

      otlt_storage(line_no)%interior_face_centre = otlt_storage(line_no)%centre &
                                                   - otlt_storage(line_no)%vein_length*otlt_storage(line_no)%outward_unit_normal

    end do

    close (unit=92)

  end subroutine setup_otlt_storage

  subroutine delete_otlt_storage()

    implicit none

    deallocate (otlt_storage)

    if (allocated(element_outlet_no)) then
      element_outlet_no = 0
      deallocate (element_outlet_no)
    end if

  end subroutine delete_otlt_storage

  subroutine setup_basis_map()

    use param

    implicit none

    integer :: cavity_no
    real(db) :: n_x, n_y, n_z, denom

    do cavity_no = 1, no_cavities

      n_x = cvty_storage(cavity_no)%orientation_normal(1)
      n_y = cvty_storage(cavity_no)%orientation_normal(2)
      n_z = cvty_storage(cavity_no)%orientation_normal(3)

      denom = sqrt(1.0_db + (n_x/n_z)**2)

      cvty_storage(cavity_no)%basis_map(1, 1) = 1.0_db/denom
      cvty_storage(cavity_no)%basis_map(1, 2) = 0.0_db
      cvty_storage(cavity_no)%basis_map(1, 3) = -(n_x/n_z)/denom
      cvty_storage(cavity_no)%basis_map(2, 1) = -(n_y*n_x/n_z)/denom
      cvty_storage(cavity_no)%basis_map(2, 2) = n_z*denom
      cvty_storage(cavity_no)%basis_map(2, 3) = -n_y/denom
      cvty_storage(cavity_no)%basis_map(3, 1) = n_x
      cvty_storage(cavity_no)%basis_map(3, 2) = n_y
      cvty_storage(cavity_no)%basis_map(3, 3) = n_z

    end do

  end subroutine setup_basis_map

  function map_global_point_to_local_cavity_coords(problem_dim, cavity_no, point)

    implicit none

    integer, intent(in) :: problem_dim, cavity_no
    real(db), dimension(problem_dim), intent(in) :: point
    real(db), dimension(problem_dim) :: &
      map_global_point_to_local_cavity_coords

    real(db), dimension(problem_dim) :: origin_point

    origin_point = cvty_storage(cavity_no)%centre

    map_global_point_to_local_cavity_coords = &
      MATMUL(cvty_storage(cavity_no)%basis_map, point - origin_point)

  end function map_global_point_to_local_cavity_coords

  function calc_dist_to_cavity_surface(problem_dim, cavity_no, point)

    implicit none

    integer, intent(in) :: problem_dim, cavity_no
    real(db), dimension(problem_dim), intent(in) :: point
    real(db) :: calc_dist_to_cavity_surface

    real(db) :: mapped_point_factor
    real(db), dimension(problem_dim) :: mapped_point

    mapped_point = &
      map_global_point_to_local_cavity_coords(problem_dim, cavity_no, point)

    mapped_point_factor = sqrt( &
                          (mapped_point(1)/cvty_storage(cavity_no)%minor_axis)**2 &
                          + (mapped_point(2)/cvty_storage(cavity_no)%minor_axis)**2 &
                          + (mapped_point(3)/cvty_storage(cavity_no)%major_axis)**2)

    calc_dist_to_cavity_surface = &
      (1.0_db - 1.0_db/mapped_point_factor)*norm2(mapped_point)

  end function calc_dist_to_cavity_surface

  function cavity_ellipsoid_eval(problem_dim, cavity_no, point)

    implicit none

    integer, intent(in) :: problem_dim, cavity_no
    real(db), dimension(problem_dim), intent(in) :: point
    real(db) :: cavity_ellipsoid_eval

    ! Local vars
    real(db) :: minor_axis, major_axis, minor_component, major_component
    real(db), dimension(problem_dim) :: mapped_point, origin_point

    minor_axis = cvty_storage(cavity_no)%minor_axis
    major_axis = cvty_storage(cavity_no)%major_axis

    mapped_point = map_global_point_to_local_cavity_coords(problem_dim, cavity_no, point)

    minor_component = (mapped_point(1)**2 + mapped_point(2)**2)/minor_axis**2
    major_component = (mapped_point(3)**2)/major_axis**2

    ! Essentially a curvilinear distance function from surface of ellipsoid
    cavity_ellipsoid_eval = minor_component + major_component - 1.0_db

  end function cavity_ellipsoid_eval

  integer function bdry_no_to_bdry_index(bdry_no)

    implicit none

    integer, intent(in) :: bdry_no

    bdry_no_to_bdry_index = bdry_no - 1

  end function bdry_no_to_bdry_index

end module boundary_data_storage
