module coupled_data_storage

  use coupled_data_type

  integer, parameter :: coupled_solns = 1

  type(fe_data_structure), dimension(coupled_solns) :: fe_data_struct

  type(coupled_user_data) :: scheme_coupled_data

end module coupled_data_storage
