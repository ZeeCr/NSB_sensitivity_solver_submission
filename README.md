# NSB_sensitivity (SWIRL)

This repository contains Fortran codes used for solving the flow and oxygen transport within placental geometries as part of the Wellcome Leap In Utero (SWIRL) project. The PDEs solved for this are unsteady, semi-discretised (BDF) Navier-Stokes-Brinkman (NSB) and steady advection-diffusion-reaction (ADR), respectively. The PDEs and associated biomarker calculations are separated into distinct programs within the subdirectories. Each subdirectory contains Fortran code and a Makefile (running `make` is sufficient within each subdirectory) for building the corresponding executable.

## Summary

- Language: Fortran (2003 standard)
- Purpose: numerical solution of flow and steady ADR for blood and oxygen transport within placental geometries and calculation of relevant biomarkers
- Build: each subdirectory contains a `Makefile` to build its executable (run `make` in each subdirectory)

## Repository layout

- `flow/` - solves NSB in (fixed) placental geometries
  - `main.f90` - main program for NSB solver
  - `aptofem_control_file.dat` - user-defined and FE parameter initialisation file
  - `ns_jac_matrix_residual_bdf_nondim.f90` - NSB assembly
  - `Makefile` - simple 'make' setup
  - `coupled_data_storage.f90` - declaration of coupled FE objects

- `markers/` - calculates flow-relevant markers
  - `main.f90` - main program with marker calculations

- `ADR_steady/` - solves steady ADR in placental geometries using the solution of `flow/`
  - `adr_av_jac_matrix_residual.f90` - ADR assembly

- `ADR_steady_markers/`

- Top-level files:
  - `flow_ADR_steady.py` - Python automation script, strongly recommended to run code through this as it handles automatic setting of certain parameters and file names
  - `startup.py` - screen and log setup script
  - `README.md` - README
  - `NSB_guide.pdf` - A PDF with a summary of steps to run this code (assumes prerequisites are already installed)

## Requirements

- gfortran >= 4.9 (or equivalent Fortran compiler with support for 2003 standard)
- AptoFEM finite element library (for further information or access, contact Paul.Houston@nottingham.ac.uk)
- an existing `data/` directory in the root directory corresponding to a particular idealised placental geometry as created by the Gmsh/OCC/Python idealised placenta code

### Optional

- (Recommended) Open MPI >= 4.1.0

## Build and run

This project uses simple Makefiles located in each subdirectory. Typical build steps (run inside a subdirectory):

```pwsh
cd {subdir}
make
(mpirun -n X) ./{program}.out
```

Check the `output/` directory for VTK solution files.

## Notes

The Fortran code is hard-coded to read in files with specific names. In particular:
- `main.f90` of flow marker, it expects a cardiac cycle to be covered by 70 timesteps.
- flow marker, ADR steady and ADR markers expect the aptofem_run_no from the relevant calculations (therefore saved solution files) to be 2.

The cardiac cycle requires small changes to faciliate other timesteps. The aptofem_run_no is automatically handled in the Python script, otherwise manual changes to the hard-coded values on solution read are needed.

## Example solution fields

The colour bars below correspond to flow speed. Images have been resized for readability.

### Local (central cavity) recirculating flow visualisation

<img src="https://github.com/user-attachments/assets/829aa146-ef8d-4bc0-8b4f-1ede7e014ec1" alt="Local recirculating flow" width="480" />

### Top-down view of placental geometry with path lines via particle tracking

<img src="https://github.com/user-attachments/assets/48ae1125-39cc-45e4-9ff1-70268b9df0e6" alt="Top-down view with path lines" width="480" />

### As above with side-on view

<img src="https://github.com/user-attachments/assets/71ff0b4e-6cb4-4424-a67e-f0a67b2ac1ae" alt="Side-on view" width="480" />

### Flow field in placental geometry with a non-symmetrical moving boundary (contraction)

<img src="https://github.com/user-attachments/assets/b98cd791-1ed4-43a5-9a0e-0212b244475d" alt="Non-symmetrical moving boundary" width="480" />