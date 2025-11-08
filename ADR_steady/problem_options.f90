module problem_options

  use aptofem_kernel

  save

! User data

  real(db) :: interior_penalty_parameter = 10.0_db !< Size of the
  !< integer penalty parameter (DGFEM only)
  character(len=aptofem_length_key_def) :: solver_section, fe_solution_section, &
    adv_solution_section

  real(db) :: Pe
  real(db) :: Dm

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

    call get_aptofem_key_definition('solver_section', &
                                    solver_section, section_name, aptofem_stored_keys, ierr)

    if (ierr /= 0) then
      write (io_err, *) 'Error: get_user_data'
      write (io_err, *) 'solver_section is not specified'
      stop
    end if

    call get_aptofem_key_definition('fe_solution_section', &
                                    fe_solution_section, section_name, aptofem_stored_keys, ierr)

    if (ierr /= 0) then
      write (io_err, *) 'Error: get_user_data'
      write (io_err, *) 'fe_solution_section is not specified'
      stop
    end if

    call get_aptofem_key_definition('adv_solution_section', &
                                    adv_solution_section, section_name, aptofem_stored_keys, ierr)

    if (ierr /= 0) then
      write (io_err, *) 'Error: get_user_data'
      write (io_err, *) 'adv_solution_section is not specified'
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
      write (io_warn, *) 'Error: get_user_data'
      write (io_warn, *) 'Dm is not specified'
      write (io_warn, *) 'Setting Dm = 0'
      Dm = 0.0_db
    end if

  end subroutine get_user_data
!--------------------------------------------------------------------
  integer function replace_key_val_acf_real_db(cf_name, cf_key, cf_val)

    use param
    use aptofem_kernel

    implicit none

    character(len=*), intent(in) :: cf_name, cf_key
    real(db), intent(in) :: cf_val
    character(len=len(cf_name) + 5) :: cf_temp_name
    character(len=256) :: line, altered_line, val_str
    character(len=8) :: advance_flag
    integer :: unit, ierr, line_no, no_lines, no_entries

    no_lines = 0
    no_entries = 0

    if (processor_no == 0) then
      unit = 201
      open (unit, file=trim(adjustl(cf_name)), status='old', action='read')
      do
        read (unit, '(A)', iostat=ierr) line
        if (ierr == 0) then
          no_lines = no_lines + 1
        else
          exit
        end if
        if (index(line, trim(adjustl(cf_key))) == 1) then
          no_entries = no_entries + 1
        end if
      end do

      ! Error check
      if (no_entries /= 1) then
        call write_message(io_err, 'Error reading '//trim(adjustl(cf_name)))
        call write_message(io_err, 'Number keys found ', no_entries)
        replace_key_val_acf_real_db = no_entries
        STOP no_entries
      end if

      ! Fortran decides in sequential files to remove every line past existing line when writing
      ! So need to create new file
      cf_temp_name = 'temp_'//trim(adjustl(cf_name))
      open (unit + 1, file=trim(adjustl(cf_temp_name)), action='write')
      write (val_str, '(D16.8)') cf_val
      write (altered_line, '(A)') trim(adjustl(cf_key))//' '//trim(adjustl(val_str))
      rewind (unit)
      do line_no = 1, no_lines
        read (unit, '(A)', iostat=ierr) line
        if (ierr /= 0) exit
        if (line_no == no_lines) then
          advance_flag = 'no'
        else
          advance_flag = 'yes'
        end if
        if (index(line, trim(adjustl(cf_key))) == 1) then
          write (unit + 1, '(A)', advance=trim(adjustl(advance_flag))) trim(adjustl(trim(altered_line)))
        else
          write (unit + 1, '(A)', advance=trim(adjustl(advance_flag))) trim(adjustl(line))
        end if
      end do
      close (unit)
      close (unit + 1)

      ierr = remove_file(cf_name)
      if (ierr /= 0) then
        call write_message(io_err, 'Error removing file: ', cf_name)
        STOP - 1
      end if

      ! Rename the new file
      ierr = rename_file(cf_temp_name, cf_name)
      if (ierr /= 0) then
        call write_message(io_err, 'Error renaming file: ', cf_temp_name)
        STOP - 1
      end if

      replace_key_val_acf_real_db = 0

    end if

#ifdef MPI
    call MPI_Barrier(mpi_communicator, ierr)
    call MPI_Bcast(replace_key_val_acf_real_db, 1, MPI_INTEGER, 0, mpi_communicator, ierr)
#endif

  end function replace_key_val_acf_real_db

  integer function remove_file(file)

    implicit none

    character(len=*), intent(in) :: file
    character(len=256) :: command

    command = 'rm '//trim(adjustl(file))

    remove_file = system(command)

  end function remove_file

  integer function rename_file(old_file, new_file)

    implicit none

    character(len=*), intent(in) :: old_file, new_file
    character(len=256) :: command

    command = 'mv '//trim(adjustl(old_file))//' '//trim(adjustl(new_file))

    rename_file = system(command)

  end function rename_file

end module problem_options
