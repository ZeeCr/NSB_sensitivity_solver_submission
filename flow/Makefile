include $(AptoFEM_DIR)/source/package_defaults/makefile_mpi.flags

main_sources = problem_options.f90 \
	input_data_storage.f90 boundary_data_storage.f90 coupled_data_storage.f90 \
	geometry_routines.f90 \
    general_routines.f90 \
    bdry_interpolate.f90 bdry_map.f90 \
	ns_jac_matrix_residual_bdf_nondim.f90

ns_bdf ::
	($(FC) $(FFLAGS) $(MODFLAGS) -o $@.out \
        $(main_sources) bdry_file.f90 main.f90 \
        $(AptoFEM) -lpthread)

cleanall ::
	(rm -f *.o)
	(rm -f *~)
	(rm -f *.mod)
	(rm -f *.out)

