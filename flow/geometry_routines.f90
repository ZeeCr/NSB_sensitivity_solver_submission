module geometry_routines

  use param
  use linear_algebra

  use boundary_data_storage

contains

  subroutine set_element_cavity_no(mesh_data)

    use fe_mesh
    use functions

    type(mesh), intent(in) :: mesh_data

    ! Local vars
    integer :: problem_dim, no_eles, ele_no, cavity_no, outlet_no
    real(db) :: minor_axis, major_axis, minor_component, major_component, &
                surface_dist, centre_dist, min_dist, x, y, z, dist_to_line, mapped_point_factor
    real(db), dimension(mesh_data%problem_dim) :: centroid, ptA, ptB, AB, AxB, &
                                                  closest_pt_on_line, perp_line, mapped_point
    real(db), dimension(2) :: bary_coords, coord_rhs
    real(db), dimension(2, 2) :: coord_mat

    real(db) :: tol = 1.0e-6

    problem_dim = mesh_data%problem_dim
    no_eles = mesh_data%no_eles

    allocate (element_cavity_no(no_eles), element_outlet_no(no_eles))
    element_cavity_no = 0
    element_outlet_no = 0

    do ele_no = 1, no_eles

      call compute_ele_centroid(centroid, problem_dim, &
                                ele_no, mesh_data)

      x = centroid(1)
      y = centroid(2)
      z = centroid(3)

      !!!! HARD-CODED - RADIUS OF GEOMETRY BUT I DON'T STORE THIS YET, THE VALUE BELOW 14.518876558339993_db

      ! This checks whether the point is outside the main volume of the placenta, i.e. in inlets/outlets
      if (x**2 + y**2 + (z - 14.518876558339993_db)**2 &
          > (14.518876558339993_db - tol)**2) then

        element_cavity_no(ele_no) = -1
        element_outlet_no(ele_no) = -1

        ! More in-depth search
      else

        ! Checks whether inside a cavity: if surface_dist < 0, then yes inside one
        do cavity_no = 1, no_cavities

          surface_dist = cavity_ellipsoid_eval(problem_dim, cavity_no, centroid)
          if (surface_dist < tol) then
            element_cavity_no(ele_no) = -1
            element_outlet_no(ele_no) = -1
            exit
          end if

        end do

      end if

      ! So point is somewhere inside bulk of placenta, including septal wall veins
      if (element_cavity_no(ele_no) == 0) then

        ! First check whether in any of the outlets
        do outlet_no = 1, no_outlets

          ! A = outer face point, B = point on placenta wall
          ptA = otlt_storage(outlet_no)%centre
          ptB = otlt_storage(outlet_no)%interior_face_centre

          AB = ptB - ptA
          call vec_product(ptA, ptB, AxB)
          call vec_product(AB, centroid, perp_line)
          dist_to_line = norm2(AxB + perp_line)/norm2(AB)

          if (dist_to_line < otlt_storage(outlet_no)%total_radius + tol) then

            call vec_product(AB, AxB + perp_line, closest_pt_on_line)
            closest_pt_on_line = centroid + closest_pt_on_line/dot_product(AB, AB)

            coord_mat(1, 1) = 1.0_db + dot_product(ptA, ptA)
            coord_mat(1, 2) = 1.0_db + dot_product(ptA, ptB)
            coord_mat(2, 1) = coord_mat(1, 2)
            coord_mat(2, 2) = 1.0_db + dot_product(ptB, ptB)

            coord_rhs(1) = 1.0_db + dot_product(closest_pt_on_line, ptA)
            coord_rhs(2) = 1.0_db + dot_product(closest_pt_on_line, ptB)

            bary_coords = solve_ax_eq_b(coord_mat, coord_rhs, 2)

            if (bary_coords(1) > -tol .AND. bary_coords(1) < 1.0_db + tol .AND. &
                bary_coords(2) > -tol .AND. bary_coords(2) < 1.0_db + tol) then

              element_cavity_no(ele_no) = -1
              element_outlet_no(ele_no) = -1
              exit

            end if

          end if

        end do

        ! Final step - must be in bulk so see which cavity / outlet it's closest to
        ! Note called centre dist, but actually distance from cavity
        ! Also veins have an implicit cavity, it's circle with axis given by funnel radius
        if (element_cavity_no(ele_no) == 0) then

          min_dist = 1.0e10
          do cavity_no = 1, no_cavities

            centre_dist = calc_dist_to_cavity_surface(problem_dim, cavity_no, centroid)

            if (centre_dist < min_dist) then
              element_cavity_no(ele_no) = cavity_no
              min_dist = centre_dist
              element_outlet_no(ele_no) = -2
            end if

          end do
          do outlet_no = 1, no_outlets

            centre_dist = norm2(centroid - otlt_storage(outlet_no)%interior_face_centre)

            if (centre_dist < min_dist) then
              element_outlet_no(ele_no) = outlet_no
              min_dist = centre_dist
              element_cavity_no(ele_no) = -2
            end if

          end do

        end if

      end if

    end do

    ! Final check
    do ele_no = 1, no_eles
      if (element_cavity_no(ele_no) == 0 .OR. element_outlet_no(ele_no) == 0 .OR. &
          (element_cavity_no(ele_no) == -2 .AND. element_outlet_no(ele_no) == -2)) then
        print *, "ERROR: set_element_cavity_no"
        print *, "Element has value 0 or both -2"
        STOP - 1
      end if
    end do

  end subroutine set_element_cavity_no

!--------------------------------------------------------------------
  subroutine eval_tanh_smooth( &
    problem_dim, ele_no, point, perm, eff_perm)
!--------------------------------------------------------------------

    implicit none

    integer, intent(in) :: problem_dim, ele_no
    real(db), intent(in) :: perm !< Exterior permeability
    real(db), dimension(problem_dim) :: point
    real(db), intent(out) :: eff_perm

    integer :: outlet_no
    real(db) :: dist

    if (element_cavity_no(ele_no) == -1 .OR. element_outlet_no(ele_no) == -1) then
      dist = 0.0_db
    else if (element_cavity_no(ele_no) == -2) then
      outlet_no = element_outlet_no(ele_no)
      dist = norm2(point - otlt_storage(outlet_no)%interior_face_centre)
      ! HACK TO CREATE A DOME FOR NOW
      if (dist > otlt_storage(outlet_no)%total_radius) then
        dist = dist - otlt_storage(outlet_no)%total_radius
      else
        dist = 0.0_db
      end if
    else if (element_outlet_no(ele_no) == -2) then
      ! In general_routines, set_element_cavity_no() already sets up the array
      ! so that element_cavity_no(ele_no) is the cavity no. of closest cavity
      ! to element ele_no
      !dist = cavity_ellipsoid_eval(problem_dim,element_cavity_no(ele_no),point)
      dist = &
        calc_dist_to_cavity_surface(problem_dim, element_cavity_no(ele_no), point)
    else
      call write_message(io_err, "ERROR: calc_eff_permeability_reciprocal")
      call write_message(io_err, "ERROR: case I didn't think of")
      write (io_err, *) element_cavity_no(ele_no), element_outlet_no(ele_no)
    end if

    eff_perm = perm*tanh(smoothing_freq*dist)

  end subroutine eval_tanh_smooth

  subroutine read_no_outlets(no_basal, no_peripheral, no_septal)

    use param

    implicit none

    character(len=256) line
    integer :: io_code, line_no, ipos1, ipos2, counter
    integer, intent(out) :: no_basal, no_peripheral, no_septal

    no_basal = 0
    no_peripheral = 0
    no_septal = 0

    open (unit=90, file='../data/geom_info.csv', iostat=io_code, status='old')

    if (io_code /= 0) then
      write (io_err, *) 'Reading file geom_info.csv failed: ', io_code
      stop
    end if

    ! Skip header
    read (90, *)

    read (90, '(A256)') line

    ipos1 = 1
    ipos2 = 2
    counter = 0
    do while (ipos2 > 1)
      ipos2 = index(line(ipos1:), ',')

      if (counter == 9) then
        read (line(ipos1:ipos1 + ipos2 - 1), *) no_basal
      else if (counter == 10) then
        read (line(ipos1:ipos1 + ipos2 - 1), *) no_peripheral
      else if (counter == 11) then
        read (line(ipos1:), *) no_septal
      end if

      ipos1 = ipos1 + ipos2
      counter = counter + 1
    end do

    call write_message(io_msg, "no_basal "//num2str(no_basal))
    call write_message(io_msg, "no_peripheral "//num2str(no_peripheral))
    call write_message(io_msg, "no_septal "//num2str(no_septal))

    close (unit=90)

  end subroutine read_no_outlets

  subroutine read_outlet_radii(no_basal, no_peripheral, no_septal, &
                               basal_r, peripheral_r, septal_r)

    use param

    implicit none

    character(len=256) line
    integer :: io_code, line_no, ipos1, ipos2, counter
    integer, intent(in) :: no_basal, no_peripheral, no_septal
    real(db), dimension(no_basal), intent(out) :: basal_r
    real(db), dimension(no_peripheral), intent(out) :: peripheral_r
    real(db), dimension(no_septal), intent(out) :: septal_r

    open (unit=90, file='../data/outlet_info.csv', iostat=io_code, status='old')

    if (io_code /= 0) then
      write (io_err, *) 'Reading file outlet_info.csv failed: ', io_code
      stop
    end if

    ! Skip header
    read (90, *)

    do line_no = 1, no_basal
      read (90, '(A256)') line

      ipos1 = 1
      ipos2 = 2
      counter = 0
      do while (ipos2 > 1)
        ipos2 = index(line(ipos1:), ',')

        if (counter == 1) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) basal_r(line_no)
          exit
        end if

        ipos1 = ipos1 + ipos2
        counter = counter + 1
      end do
    end do

    do line_no = 1, no_peripheral
      read (90, '(A256)') line

      ipos1 = 1
      ipos2 = 2
      counter = 0
      do while (ipos2 > 1)
        ipos2 = index(line(ipos1:), ',')

        if (counter == 1) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) peripheral_r(line_no)
          exit
        end if

        ipos1 = ipos1 + ipos2
        counter = counter + 1
      end do
    end do

    do line_no = 1, no_septal
      read (90, '(A256)') line

      ipos1 = 1
      ipos2 = 2
      counter = 0
      do while (ipos2 > 1)
        ipos2 = index(line(ipos1:), ',')

        if (counter == 1) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) septal_r(line_no)
          exit
        end if

        ipos1 = ipos1 + ipos2
        counter = counter + 1
      end do
    end do

    print *, "basal_r", basal_r
    print *, "peripheral_r", peripheral_r
    print *, "septal_r", septal_r

    close (unit=90)

  end subroutine read_outlet_radii

  subroutine read_no_arteries(no_arteries)

    use param

    implicit none

    character(len=256) line
    integer :: io_code, line_no, ipos1, ipos2, counter
    integer, intent(out) :: no_arteries

    no_arteries = 0

    open (unit=90, file='../data/mesh_info.csv', iostat=io_code, status='old')

    if (io_code /= 0) then
      write (io_err, *) 'Reading file mesh_info.csv failed: ', io_code
      stop
    end if

    ! Skip header
    read (90, *)

    do
      read (90, '(A256)', iostat=io_code) line
      if (io_code /= 0) then
        exit
      else
        no_arteries = no_arteries + 1
      end if
    end do

    call write_message(io_msg, "no_arteries "//num2str(no_arteries))

    close (unit=90)

  end subroutine read_no_arteries

  subroutine read_artery_radii(no_arteries, artery_r)

    use param

    implicit none

    character(len=256) line
    integer :: io_code, line_no, ipos1, ipos2, counter
    integer, intent(in) :: no_arteries
    real(db), dimension(no_arteries), intent(out) :: artery_r

    open (unit=90, file='../data/mesh_info.csv', iostat=io_code, status='old')

    if (io_code /= 0) then
      write (io_err, *) 'Reading file mesh_info.csv failed: ', io_code
      stop
    end if

    ! Skip header
    read (90, *)

    do line_no = 1, no_arteries
      read (90, '(A256)') line

      ipos1 = 1
      ipos2 = 2
      counter = 0
      do while (ipos2 > 1)
        ipos2 = index(line(ipos1:), ',')

        if (counter == 1) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) artery_r(line_no)
          exit
        end if

        ipos1 = ipos1 + ipos2
        counter = counter + 1
      end do
    end do

    print *, "artery_r", artery_r

    close (unit=90)

  end subroutine read_artery_radii

  ! Very unoptimised, but only called once
  ! This assumes that the cell nodes are ordered in a counter/clockwise order in the file
  ! Is a problem_dim array but the problem_dim index is not read, it's just 0
  subroutine read_cotyledon_node_locations(problem_dim, &
                                           no_cells, nodes_per_cell, max_nodes_per_cell, c_nodes)

    implicit none

    integer, intent(in) :: problem_dim
    integer, intent(out) :: no_cells, max_nodes_per_cell
    integer, dimension(:), allocatable, intent(inout) :: nodes_per_cell
    real(db), dimension(:, :, :), allocatable, intent(inout) :: c_nodes

    character(len=256) line
    integer :: io_code, total_nodes, line_no, ipos1, ipos2, counter, temp_int, c_no

    no_cells = 0
    max_nodes_per_cell = 0

    open (unit=90, file='../data/c_loc.csv', iostat=io_code, status='old')

    if (io_code /= 0) then
      write (io_err, *) 'Reading file c_loc.csv failed: ', io_code
      stop - 1
    end if

    read (90, *)

    ! Read number of nodes and number of cells
    total_nodes = 0
    do
      ipos1 = 1
      ipos2 = 2
      read (90, '(A256)', iostat=io_code) line
      if (io_code /= 0) then
        exit
      else
        total_nodes = total_nodes + 1
        ipos2 = index(line(ipos1:), ',')
        read (line(ipos1:ipos1 + ipos2 - 1), *) temp_int
        no_cells = max(no_cells, temp_int)
      end if
    end do

    ! Allocate number of nodes storage
    allocate (nodes_per_cell(no_cells))
    nodes_per_cell = 0

    ! Read number of nodes per cell
    rewind (unit=90)
    read (90, *) ! Skip header
    do line_no = 1, total_nodes
      ipos1 = 1
      read (90, '(A256)') line
      ipos2 = index(line(ipos1:), ',')
      read (line(ipos1:ipos1 + ipos2 - 1), *) temp_int
      nodes_per_cell(temp_int) = nodes_per_cell(temp_int) + 1
    end do
    max_nodes_per_cell = maxval(nodes_per_cell)

    ! Allocate cotyledon node storage
    allocate (c_nodes(problem_dim, max_nodes_per_cell, no_cells))
    c_nodes = 0.0_db

    rewind (unit=90)
    read (90, *) ! Skip header
    do line_no = 1, total_nodes
      read (90, '(A256)') line

      ipos1 = 1
      ipos2 = 2
      counter = 0
      do while (ipos2 > 1)
        ipos2 = index(line(ipos1:), ',')

        ! Cell number
        if (counter == 0) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) c_no
          counter = counter + 1
          ! Edge / node number
        else if (counter == 1) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) temp_int
          counter = counter + 1
          ! x
        else if (counter == 2) then
          read (line(ipos1:ipos1 + ipos2 - 1), *) c_nodes(problem_dim - 2, temp_int, c_no)
          counter = counter + 1
          ! y
        else if (ipos2 == 0 .AND. counter == 3) then
          read (line(ipos1:), *) c_nodes(problem_dim - 1, temp_int, c_no)
        else
          write (io_err, *) "Error: read_cotyledon_node_locations"
          stop - 10
        end if

        ipos1 = ipos1 + ipos2

      end do

      ! z
      c_nodes(problem_dim, temp_int, c_no) = 0.0_db

    end do

    close (unit=90)

  end subroutine read_cotyledon_node_locations

end module geometry_routines
