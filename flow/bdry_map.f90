module boundary_mapping

  implicit none

contains

!--------------------------------------------------------------------
! Allocates and maps the standard array indexing to the boundary faces
!--------------------------------------------------------------------
  subroutine boundary_map(no_bdry_map, bdry_map)

    use param
    implicit none

    integer, intent(in) :: no_bdry_map !< Number of boundaries to map to
    !< standard array indexing
    integer, dimension(no_bdry_map), intent(inout) :: bdry_map

    integer :: i

    do i = 1, no_bdry_map
      bdry_map(i) = i + 1
    end do

  end subroutine boundary_map

!--------------------------------------------------------------------
! Maps boundary number back to boundary_map array
!--------------------------------------------------------------------
  function inverse_boundary_map(bdry_no)

    use param
    implicit none

    integer, intent(in) :: bdry_no
    integer :: inverse_boundary_map

    write (io_err, *) "ERROR: inverse_boundary_map"
    write (io_err, *) "Should not be called here"
    stop - 1

    select case (bdry_no)
    case (2)
      inverse_boundary_map = 1
    case (3)
      inverse_boundary_map = 2
    case default
      write (io_err, *) "Incorrect bdry_no in inverse_boundary_map"
      write (io_err, *) "bdry_no: ", bdry_no
      stop
    end select

  end function inverse_boundary_map

end module boundary_mapping
