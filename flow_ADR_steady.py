import sys
import os
import subprocess

def main():
  
  # Linear scaling on the nominal char. vel.
  no_vel_variations = 2
  
  nominal_cardiac_cycle = 10.5 # 0.7s / (0.01 / 0.15) == cardiac_cycle (s) / (L_c / V_c)

  DaRecipr_arr = [10.0**i for i in range(4,5)]
  Re_arr = [454.286*(i+1) for i in range(0,no_vel_variations)]
  Pe_arr = [9.0e5*(i+1) for i in range(0,no_vel_variations)]
  T_arr = [float(i+1) for i in range(0,no_vel_variations)]
  
  for vel_counter in range(0,no_vel_variations):
    
    Re = Re_arr[vel_counter]
    Pe = Pe_arr[vel_counter]
    T_scaling = T_arr[vel_counter]
    
    # Only vary V_c, L_c constant
    cardiac_cycle_t = nominal_cardiac_cycle * T_scaling
    total_time = 2.0 * cardiac_cycle_t
    total_no_timesteps = 140
    time_step = total_time / total_no_timesteps
        
    subprocess.run(["date"])
    subprocess.run(["echo", f"Solving Re = {Re} and Pe = {Pe}"])
        
    for DaRecipr in DaRecipr_arr:

      subprocess.run(["date"])
      subprocess.run(["echo", f"Solving Da_reciprocal = {DaRecipr}"])

      # Read CSV file
      apto_run_no = os.path.join('flow', 'aptofem_run_number.dat')
      with open(apto_run_no, 'w') as file:
        file.write("1")
        file.close()
    
      # Update line in acf.dat file
      acf_file_path = os.path.join('flow', 'aptofem_control_file.dat')
      with open(acf_file_path, 'r') as file:
        lines = file.readlines()
        file.close()
          
      lines = replace_str_in_lines('Re', Re, lines)
      lines = replace_str_in_lines('Da_reciprocal', DaRecipr, lines)
      lines = replace_str_in_lines('bdf_final_time', total_time, lines)
      lines = replace_str_in_lines('t_scaling', T_scaling, lines)
      
      with open(acf_file_path, 'w') as file:
        file.writelines(lines)
        file.close()
          
      # Update line in acf.dat file
      acf_file_path = os.path.join('markers', 'aptofem_control_file.dat')
      with open(acf_file_path, 'r') as file:
        lines = file.readlines()
        file.close()
          
      lines = replace_str_in_lines('Re', Re, lines)
      lines = replace_str_in_lines('Da_reciprocal', DaRecipr, lines)
      lines = replace_str_in_lines('time_step', time_step, lines)
      
      with open(acf_file_path, 'w') as file:
        file.writelines(lines)
        file.close()
      
      # Run bash command
      subprocess.run(["mpirun -n 4 ./ns_bdf.out"], shell=True, cwd="./flow/")
      subprocess.run(["./markers.out"], shell=True, cwd="./markers/")
      
      subprocess.run([f"cp ./meshandsoln_soln_1_navier_2_81.vtk ./Re{Re}_DaRecipr{DaRecipr}_81.vtk"], shell=True, cwd="./flow/output/")
      subprocess.run([f"cp ./solution_for_restart_soln_1_2_81.internal ./Re{Re}_DaRecipr{DaRecipr}_81.internal"], shell=True, cwd="./flow/output/restart/")
      
      ADR(Re,Pe,DaRecipr)
    
  return None

def ADR(Re,Pe,DaRecipr):

    counter = 1
    multiple = 1.11e-5
    while (multiple < 0.112):
      
        Dm = multiple*counter
        
        subprocess.run(["date"])
        subprocess.run(["echo", f"Solving Dm = {Dm} with Re = {Re}, Pe = {Pe}, Da_recripcal = {DaRecipr}"])

        # Read CSV file
        apto_run_no = os.path.join('ADR_steady', 'aptofem_run_number.dat')
        with open(apto_run_no, 'w') as file:
            file.write("1")
            file.close()

        # Update line in acf.dat file
        acf_file_path = os.path.join('ADR_steady', 'aptofem_control_file.dat')
        with open(acf_file_path, 'r') as file:
            lines = file.readlines()
            file.close()
            
        lines = replace_str_in_lines('Pe', Pe, lines)
        lines = replace_str_in_lines('Dm', Dm, lines)
        
        with open(acf_file_path, 'w') as file:
            file.writelines(lines)
            file.close()
            
        # Update line in acf.dat file
        acf_file_path = os.path.join('ADR_steady_markers', 'aptofem_control_file.dat')
        with open(acf_file_path, 'r') as file:
            lines = file.readlines()
            file.close()
            
        lines = replace_str_in_lines('Re', Re, lines)
        lines = replace_str_in_lines('Pe', Pe, lines)
        lines = replace_str_in_lines('Dm', Dm, lines)
        lines = replace_str_in_lines('Da_reciprocal', DaRecipr, lines)
        
        with open(acf_file_path, 'w') as file:
            file.writelines(lines)
            file.close()
        
        # Run bash command
        subprocess.run(["mpirun -n 4 ./adr_steady.out"], shell=True, cwd="./ADR_steady/")
        subprocess.run(["./adr_steady_markers.out"], shell=True, cwd="./ADR_steady_markers/")
        
        if (abs(Dm - 1.11e-3) <= 1.0e-7):
            subprocess.run([f"cp ./meshandsoln_ADR_2_1.vtk ./Re{Re}_Pe{Pe}_DaRecipr{DaRecipr}.vtk"], shell=True, cwd="./ADR_steady/output/")
            
        counter = counter + 5
        if (counter > 9):
            counter = 1
            multiple = multiple*10.0
        
    return None

def replace_str_in_lines(str2find, val2replace, lines):
    line_index = None
    for i, line in enumerate(lines):
        if line.startswith(str2find):
            line_index = i
            break
    if line_index is not None:
        lines[line_index] = f'{str2find} {val2replace}\n'
    else:
        subprocess.run(["echo", f"Error: could not find {str2find}"])
        print(f"Error: could not find {str2find}")
        sys.exit(-2)
        
    return lines

def make_probs():
    
  subprocess.run(["make"], shell=True, cwd="./flow/")
  subprocess.run(["make"], shell=True, cwd="./markers/")
  subprocess.run(["make"], shell=True, cwd="./ADR_steady/")
  subprocess.run(["make"], shell=True, cwd="./ADR_steady_markers/")

try:
    make_probs()
    main()
except KeyboardInterrupt:
    subprocess.run(["echo", "Keyboard interrupt detected. Exiting..."])
    sys.exit(0)