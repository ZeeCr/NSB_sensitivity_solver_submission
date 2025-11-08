module problem_options

  use aptofem_kernel

  save

! User data

  real(db) :: interior_penalty_parameter = 10.0_db !< Size of the
  !< integer penalty parameter (DGFEM only)
  character(len=aptofem_length_key_def) :: ns_solver_section
  character(len=aptofem_length_key_def) :: fe_solution_section

  real(db) :: Re, t_scaling
  real(db) :: Da_reciprocal = 0.0_db
  real(db) :: backflow_stab_parameter = 0.0_db

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

    call get_aptofem_key_definition('Re', &
                                    re, section_name, aptofem_stored_keys, ierr)

    if (ierr /= 0) then
      write (io_err, *) 'Error: get_user_data'
      write (io_err, *) 'Re is not specified'
      stop
    end if

    call get_aptofem_key_definition('Da_reciprocal', &
                                    Da_reciprocal, section_name, aptofem_stored_keys, ierr)

    if (ierr /= 0) then
      write (io_warn, *) 'Error: get_user_data'
      write (io_warn, *) 'Da_reciprocal is not specified'
      write (io_warn, *) 'Setting Da_reciprocal = 0'
    end if

    call get_aptofem_key_definition('backflow_stab_parameter', &
                                    backflow_stab_parameter, section_name, aptofem_stored_keys, ierr)

    if (ierr /= 0) then
      write (io_warn, *) 'Warning: get_user_data'
      write (io_warn, *) 'backflow_stab_parameter not defined'
      write (io_warn, *) 'Setting backflow_stab_parameter = 0'
    end if

    call get_aptofem_key_definition('t_scaling', &
                                    t_scaling, section_name, aptofem_stored_keys, ierr)

    if (ierr /= 0) then
      write (io_warn, *) 'Warning: get_user_data'
      write (io_warn, *) 't_scaling not defined'
      write (io_warn, *) 'Setting t_scaling = 1'
      t_scaling = 1.0_db
    end if

  end subroutine get_user_data
!--------------------------------------------------------------------
end module problem_options
