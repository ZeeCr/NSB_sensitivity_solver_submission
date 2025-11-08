include $(AptoFEM_DIR)/source/package_defaults/makefile_mpi.flags

derived_sources = ../flow/boundary_data_storage.f90 \
	../flow/coupled_data_storage.f90 \
	../flow/geometry_routines.f90

main_sources = problem_options.f90 \
	adr_av_jac_matrix_residual.f90

adr_steady ::
	($(FC) $(FFLAGS) $(MODFLAGS) -o $@.out \
        $(derived_sources) $(main_sources) bdry_file.f90 main.f90 \
        $(AptoFEM) -lpthread)

cleanall ::
	(rm -f *.o)
	(rm -f *~)
	(rm -f *.mod)
	(rm -f *.out)

