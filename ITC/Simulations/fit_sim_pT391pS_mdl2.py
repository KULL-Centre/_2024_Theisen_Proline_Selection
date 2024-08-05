import numpy as np
import math
from scipy.optimize import minimize
from scipy.optimize import leastsq
from COPASI import *
import matplotlib.pyplot as plt
import datetime

# Load the model once, assuming the path to your COPASI model file
MODEL = "ITC_DataSimulation_pT391pS_mdl2"
DATAFILE = "20240612_1mM_syr_run2"
CELLVOLUME = 207.1 * 10**-6
INJVOLUME = 2.0*10**-6
INJCONC = 0
GasConstant = 0.008314
T = 310.15
N_INJ = 0

kd_cis_initial = 0.160
dg_cis_initial = -40
Hcis_initial = -20000
Htrans_initial = -20000
R_initial = 7.8
N_initial = 0.666
K_cis_trans_initial = 0.386962552
k_cis_trans_initial = 0.0004*9
K_re_initial = 0
Offset_initial = 1000

def read_data(rows):
    global N_INJ
    # Initialize a list to hold the dictionaries
    results = []
    
    # Iterate over each row
    for row in rows[1:]:
        # Split the row into parts based on commas
        parts = row.split(',')
        
        # Extract the 'cell', 'syringe', 'injmass' values (they are repeated but always the same for a given row)
        ratio = float(parts[0])
        cell = float(parts[1])
        syringe = float(parts[2])
        injmass = float(parts[3])
        v_inj = float(parts[4])
        include = parts[5] == "1"
        
        # Extract the .itc_peak values and calculate their average
        itc_peak_values = [float(parts[i]) for i in range(6, len(parts), 8)]
        avg_itc_peak = sum(itc_peak_values) / len(itc_peak_values)
        
        # Create a dictionary with the extracted values and average, and append it to the results list
        results.append({'cell': cell, 'syringe': syringe, 'injmass': injmass, 'V_inj':v_inj, 'avg_itc_peak': avg_itc_peak,'ratio':ratio, 'include':include})

    N_INJ = len(results)

    return results

def plot_model_fit(best_fit_params, actual_data):
    """
    Plots the simulated model data using best-fit parameters against the actual data.

    Parameters:
    - best_fit_params: The parameters that provide the best fit to the actual data.
    - actual_data: The actual observed data, expected to be a list of [x, y] pairs.
    """
    # Simulate the model with the best-fit parameters to get the simulated data
    simulated_data = sample_parameters(*best_fit_params)
    
    # Extracting x and y values from the actual and simulated data
    actual_x = [pair[0] for pair in actual_data]
    actual_y = [pair[1] for pair in actual_data]
    
    simulated_x = [pair[0] for pair in simulated_data]
    simulated_y = [pair[1] for pair in simulated_data]
    
    # Creating the plot
    plt.figure(figsize=(10, 6))
    plt.plot(actual_x, actual_y, 'o', label='Actual Data')
    plt.plot(simulated_x, simulated_y, 'o', label='Simulated Data', markersize=5)
    
    plt.xlabel('Molar Ratio')
    plt.ylabel('Enthalpy (kJ / mol)')
    plt.title(f'{MODEL} Fit to {DATAFILE}')
    plt.legend()
    plt.grid(True)
    plt.show()

def sample_parameters(dg_cis, Hcis, Htrans, ratio, N, K, k, Offset):
    global CELLVOLUME
    global INJCONC
    INJCONC = 0

    dataModel = CRootContainer.addDatamodel()

    # Load the model
    dataModel.loadModel('./' + MODEL + '.cps')

    # Access the model to change the parameter
    model = dataModel.getModel()
    parameters = model.getModelValues()
    kd_cis = 1000000*math.exp(dg_cis/(T*GasConstant))
    kd_trans = np.clip(ratio,1,100000) * kd_cis
    
    parameters["Kd_cis"].setInitialValue(kd_cis)
    parameters["Kd_trans"].setInitialValue(kd_trans)
    parameters["N"].setInitialValue(N)
    parameters["K_cis_trans"].setInitialValue(K)
    parameters["V_inj"].setInitialValue(INJVOLUME)
    parameters["k_cis_trans"].setInitialValue(k)
    #parameters["k_re_cis_trans"].setInitialValue(k_re)

    model.applyInitialValues() # Update values

    for m in model.getMetabolites(): # Old stuff
        name = m.getObjectName()
        if 'PRLR' in name and m.getCompartment().getObjectName() == "syringe":
            INJCONC += m.getInitialConcentration() # in umol
    
    model.applyInitialValues() # Update values

    #parameters = model.getModelValues()

    #for par in parameters:
    #    print(par,par.getValue())

    # get the trajectory task object
    trajectoryTask = dataModel.getTask("Time-Course")
    dataModel.getModel().setInitialTime(0.0)
    model.applyInitialValues() # Update values

    trajectoryTask.process(True)

    result = read_time_course_results("./output.txt")
    out = process_result(result, options = {'cis':Hcis,'trans':Htrans, 'n': N, 'offset': Offset})

    return out

# Read output file
def read_time_course_results(file_path):
    # Flag to indicate if the time-course result section has started
    parsing_results = True
    # List to store the parsed data
    results = []
    # Header for creating dictionaries
    headers = []

    # Open and read the file
    with open(file_path, 'r') as file:
        for line in file:
            # Check if the time-course result section has started
            if "Time-Course Result" in line: # Legacy code?? keep for now
                parsing_results = True
                continue
            # Check for the start of the data section within the results
            if parsing_results and line.startswith("#"):
                headers = line.strip().lstrip('#').strip().split('\t')
                continue
            # Store the data lines
            if parsing_results and line.strip():
                # Split the line into values and pair with headers
                row_data = line.strip().split('\t')
                row_dict = dict(zip(headers, row_data))
                # Convert numeric values from strings to appropriate types
                for key, value in row_dict.items():
                    try:
                        row_dict[key] = float(value)
                    except ValueError:
                        pass  # Keep as string if conversion fails
                results.append(row_dict)

    return results[1:N_INJ+2] #Skip first artifact row, second row is starting conditions so take two additional rows???

# Process output data 
def process_result(result, options):
    data = []
    prev_inj = result[0]

    for inj in result[1:]:
        INJMASS = INJCONC * prev_inj['V_inj']

        tot_prlr = sum([value * key.count('PRLR') for key, value in inj.items()]) # Can handle higher order complexes if properly named
        tot_14_3_3 = sum([value * key.count('14-3-3') for key, value in inj.items()]) / options['n']

        q_i_cis = inj['14-3-3_PRLR_Cis_Bound'] * options['cis']
        q_i_trans = inj['14-3-3_PRLR_Trans_Bound'] * options['trans']
        q_i = q_i_cis + q_i_trans

        q_i_1_cis = prev_inj['14-3-3_PRLR_Cis_Bound'] * options['cis']
        q_i_1_trans = prev_inj['14-3-3_PRLR_Trans_Bound'] * options['trans']
        q_i_1 = q_i_1_cis + q_i_1_trans

        enthalpy = q_i + (prev_inj['V_inj'] / CELLVOLUME) * ((q_i + q_i_1) / 2) - q_i_1 # kJ/mol * umolar + kJ/mol * umolar
        enthalpy /= INJMASS # moles injected, all concentrations in umol ((mJ / liter) / umol) = mJ / u * liter * mol = kJ / mol * liter
        enthalpy *= CELLVOLUME # kJ / (mol * liter) * liter = kJ / mol
        enthalpy += options['offset']

        ratio = tot_prlr/tot_14_3_3

        data.append([ratio,enthalpy])

        prev_inj = inj

    return data

def compute_rmsd(simulation_output, data_list):
    """
    Compute the RMSD between the simulation output and the provided data list.
    Both inputs are expected to be lists of [x, y] pairs.

    Parameters:
    - simulation_output: List of [x, y] pairs from the simulation.
    - data_list: List of [x, y] pairs of the reference data.

    Returns:
    - RMSD as a float.
    """
    #simulation_x_values = [pair[0] for pair in simulation_output]
    #data_x_values =  [pair[0] for pair in data_list]

    squared_differences = []

    for i in range(len(simulation_output)):
        sim = simulation_output[i]
        data = data_list[i]

        if data[2]:
            squared_differences.append((data[1] - sim[1]) ** 2)

    rmsd = np.sqrt(np.mean(squared_differences))

    return rmsd

# Define the error function for fitting
def error_function(params):
    simulated_data = sample_parameters(*params)
    error = compute_rmsd(simulated_data, DATA) # Simple sum of squared differences
    print(error)
    return error

def find_allowed_deviation(best_fit_params, actual_data, percent_increase_allowed=5):
    """
    Finds the allowed deviation for each parameter that results in an RMSD 
    up to a specified percent higher than the best fit RMSD.

    Parameters:
    - best_fit_params: Array of best-fit parameters.
    - actual_data: The actual data to compare against the simulation.
    - percent_increase_allowed: The allowed increase in RMSD, in percent.

    Returns:
    - A dictionary with parameter indices as keys and (min_deviation, max_deviation) as values.
    """
    initial_rmsd = compute_rmsd(sample_parameters(*best_fit_params), actual_data)
    allowed_rmsd = initial_rmsd * (1 + percent_increase_allowed / 100)
    
    deviations = []
    
    for i, param in enumerate(best_fit_params):
        print(f"Param: {i}, {param}")
        min_deviation = 0
        max_deviation = 0

        if param == 0: continue

        # Vary parameter upwards
        step_size = (1 + abs(param)) * 0.001  # Start with 1% of the parameter value
        new_param = param
        j = 0
        while True:
            new_param += step_size
            new_params = best_fit_params.copy()
            new_params[i] = new_param
            new_rmsd = compute_rmsd(sample_parameters(*new_params), actual_data)
            if new_rmsd > allowed_rmsd:
                break
            max_deviation = new_param - param

            j += 1
            if j % 50 == 0:
                step_size *= 2
            if j > 1000:
                print(" -no upper limit")
                break

        # Vary parameter downwards
        step_size = (1 + abs(param)) * 0.001  # Start with 1% of the parameter value
        new_param = param
        j = 0
        while True:
            new_param -= step_size
            new_params = best_fit_params.copy()
            new_params[i] = new_param
            new_rmsd = compute_rmsd(sample_parameters(*new_params), actual_data)
            if new_rmsd > allowed_rmsd:
                break
            min_deviation = new_param - param

            j += 1
            if j % 50 == 0:
                step_size *= 2
            if j > 1000:
                print(" -no lower limit")
                break

        print(f"  {param}, [{param + min_deviation}, {param + max_deviation}]")
        
        deviations.append(np.mean([abs(min_deviation), abs(max_deviation)]))

        print(f"SD = {deviations[i]}")
        print()
    
    print()
    
    return deviations

DATA = []
with open(f'{DATAFILE}.csv') as f:
    dat = read_data(f.readlines())

    for inj in dat:
        DATA.append([inj['ratio'],inj['avg_itc_peak'],inj['include']])

# Initial guesses for parameters
initial_guess = [dg_cis_initial, Hcis_initial, Htrans_initial, R_initial, N_initial, K_cis_trans_initial, k_cis_trans_initial, Offset_initial]

bounds = [
          (-70, -20),  # dG_cis should remain negative
          (-300000, 300000),  # Hcis should remain negative
          (-300000, 300000),  # Htrans has no specific bounds
          (R_initial, R_initial),      # R should be non-negative
          (N_initial*0.7, N_initial*1.5),   # N 
          (K_cis_trans_initial, K_cis_trans_initial),   # isomerization equlibrium constant
          (k_cis_trans_initial/100,k_cis_trans_initial*100),
          (-30000,30000)]

# Perform the optimization to fit the model parameters
result = minimize(error_function, initial_guess, method='Nelder-Mead', bounds=bounds, options={'adaptive':True, 'fatol': 0.1})

result = minimize(error_function, result.x, method='Nelder-Mead', bounds=bounds, options={'adaptive':False})

# Extract the fitted parameters
fitted_params = result.x
print()
print(fitted_params)
print()
print("Getting errors...")

deviations = find_allowed_deviation(fitted_params, DATA)

print()
print(f"OUTPUT: {DATAFILE}")

print(f'dG_cis = [{fitted_params[0]},{deviations[0]}]')
print(f'Hcis = [{fitted_params[1]/1000},{deviations[1]/1000}]')
print(f'Htrans = [{fitted_params[2]/1000},{deviations[2]/1000}]')
print(f'N = [{fitted_params[4]},{deviations[4]}]')
print(f'R = [{fitted_params[3]},{deviations[3]}]')
print(f'K_cis_trans = [{fitted_params[5]},{deviations[5]}]')
print(f'- k_cis_trans = [{fitted_params[6]},{deviations[6]}]')

print()

kd_cis = 10**9 * math.exp(fitted_params[0] / (T*GasConstant))
kd_cis_dist = []

kd_trans = fitted_params[3] * kd_cis
kd_trans_dist = []

dg_trans = GasConstant * T * math.log(kd_trans * 10**-9)
dg_trans_dist = []

for i in range(10000):
    _r = np.clip(np.random.normal(fitted_params[3], deviations[3], 1)[0],0.1,100000)
    _dg_cis = np.random.normal(fitted_params[0], deviations[0], 1)[0]
    _kd_cis = 10**9 * math.exp(_dg_cis / (T*GasConstant))
    kd_cis_dist.append(_kd_cis)
    kd_trans_dist.append(_r * _kd_cis)
    dg_trans_dist.append(GasConstant * T * math.log(kd_trans_dist[-1] * 10**-9))

print(f'Kd_cis = {kd_cis} ± {np.std(kd_cis_dist)} nM')
print(f'∆G_cis = {fitted_params[0]} ± {deviations[0]} kJ/mol')

print(f'Kd_trans = {kd_trans} ± {np.std(kd_trans_dist)} nM')
print(f'∆G_trans = {dg_trans} ± {np.std(dg_trans_dist)} kJ/mol')

plot_model_fit(fitted_params, DATA)

inp = input("Save Simulation?: ")

if 'y' in inp:
    out = sample_parameters(*fitted_params)
    injections = []
    for v in out:
        injections.append([v[0],[v[1]]])

    for i in range(99):
        out = sample_parameters(
            np.random.normal(fitted_params[0], deviations[0], 1)[0], 
            np.random.normal(fitted_params[1], deviations[1], 1)[0],
            np.random.normal(fitted_params[2], deviations[2], 1)[0],
            np.random.normal(fitted_params[3], deviations[3], 1)[0],
            np.clip(np.random.normal(fitted_params[4], deviations[4], 1)[0], 0.1, 10),
            np.random.normal(fitted_params[5], deviations[5], 1)[0],
            np.clip(np.random.normal(fitted_params[6], deviations[6], 1)[0],0.0001,2),
            np.random.normal(fitted_params[7],deviations[7]))

        for i in range(len(out)):
            v = out[i]
            injections[i][1].append(v[1])

    with open('output_' + MODEL + '_' + datetime.datetime.now().strftime("%d%m%Y_%H%M%S") + '.txt','w') as f:
        f.write("PARAMETERS:\n")
        f.write('dG_cis = ' + str(fitted_params[0]) + '±' + str(deviations[0]) + '\n')
        f.write('Hcis = ' + str(fitted_params[1]) + '±' + str(deviations[1]) + '\n')
        f.write('Htrans = ' + str(fitted_params[2]) + '±' + str(deviations[2]) + '\n')
        f.write('N = ' + str(fitted_params[4]) + '±' + str(deviations[4]) + '\n')
        f.write('Kd Ratio = ' + str(fitted_params[3]) + '±' + str(deviations[3]) + '\n')
        f.write('K_cis_trans = ' + str(fitted_params[5]) + '±' + str(deviations[5]) + '\n')
        f.write('k_cis_trans = ' + str(fitted_params[6]) + '±' + str(deviations[6]) + '\n')
        f.write('\n')
        f.write(f'Kd_cis = {kd_cis:.2f} ± {np.std(kd_cis_dist):.2f} nM\n')
        f.write(f'∆G_cis = {fitted_params[0]:.2f} ± {deviations[0]:.2f} kJ/mol\n')
        f.write(f'Kd_trans = {kd_trans:.2f} ± {np.std(kd_trans_dist):.2f} nM\n')
        f.write(f'∆G_trans = {dg_trans:.2f} ± {np.std(dg_trans_dist):.2f} kJ/mol\n')
        f.write('\n')
        f.write('Molar Ratio / Heats\n')
        for inj in injections:
            line = str(inj[0])

            for v in inj[1]:
                line += " " + str(v)

            f.write((line) + '\n')





