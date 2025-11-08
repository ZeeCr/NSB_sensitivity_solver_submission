module problem_options

  use aptofem_kernel

  save

! User data

  real(db) :: interior_penalty_parameter = 10.0_db !< Size of the
!< integer penalty parameter (DGFEM only)
  character(len=aptofem_length_key_def) :: ns_solver_section
  character(len=aptofem_length_key_def) :: fe_solution_section

  real(db) :: Pe, Dm, Da_reciprocal, Re

contains

!--------------------------------------------------------------------
  subroutine get_user_data(section_name, aptofem_stored_keys)
!--------------------------------------------------------------------

    implicit none

    type(aptofem_keys), pointer :: aptofem_stored_keys !< Linked
!< list used to store the AptoFEM keys
    character(len=*), intent(in) :: section_name !< Section in
!< the control file where the AptoFEM keys are defined

! Local variables

    integer :: ierr, i

    call get_aptofem_key_definition('interior_penalty_parameter', &
                                    interior_penalty_parameter, section_name, aptofem_stored_keys, ierr)

    if (ierr /= 0) then
      write (io_err, *) 'Error: get_user_data'
      write (io_err, *) 'interior_penalty_parameter is not specified'
      stop
    end if

    call get_aptofem_key_definition('ns_solver_section', &
                                    ns_solver_section, section_name, aptofem_stored_keys, ierr)

    if (ierr /= 0) then
      write (io_err, *) 'Error: get_user_data'
      write (io_err, *) 'ns_solver_section is not specified'
      stop
    end if

    call get_aptofem_key_definition('fe_solution_section', &
                                    fe_solution_section, section_name, aptofem_stored_keys, ierr)

    if (ierr /= 0) then
      write (io_err, *) 'Error: get_user_data'
      write (io_err, *) 'fe_solution_section is not specified'
      stop
    end if

    call get_aptofem_key_definition('Pe', &
                                    Pe, section_name, aptofem_stored_keys, ierr)

    if (ierr /= 0) then
      write (io_err, *) 'Error: get_user_data'
      write (io_err, *) 'Pe is not specified'
      stop
    end if

    call get_aptofem_key_definition('Dm', &
                                    Dm, section_name, aptofem_stored_keys, ierr)

    if (ierr /= 0) then
      write (io_err, *) 'Error: get_user_data'
      write (io_err, *) 'Dm is not specified'
      stop
    end if

    call get_aptofem_key_definition('Re', &
                                    Re, section_name, aptofem_stored_keys, ierr)

    if (ierr /= 0) then
      write (io_err, *) 'Warning: get_user_data'
      write (io_err, *) 'Re of underlying flow not specified'
      stop
    end if

    call get_aptofem_key_definition('Da_reciprocal', &
                                    Da_reciprocal, section_name, aptofem_stored_keys, ierr)

    if (ierr /= 0) then
      write (io_err, *) 'Error: get_user_data'
      write (io_err, *) 'Da_reciprocal is not specified'
      stop
    end if

  end subroutine get_user_data
!--------------------------------------------------------------------
end module problem_options
