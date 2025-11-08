module coupled_data_storage

  use coupled_data_type

  integer, parameter :: coupled_solns = 1

  type(fe_data_structure), dimension(coupled_solns) :: fe_data_struct

end module coupled_data_storage
