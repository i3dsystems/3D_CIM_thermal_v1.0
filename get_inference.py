# Import modules
import subprocess
import numpy as np

# Get the drift coefficient
def drift_file(var, val):
    if var==0:
        filename = 'device_binary_rram_model/output_results/results/drift_coefficient_0.txt'
    else:
        filename = 'device_binary_rram_model/output_results/results/drift_coefficient_' + str(var) + "_" + str(round(val)) + ".txt"
    with open(filename) as f:
        value = f.read()
    return value

def callNeuroSim(var, vals):
    # Call NeuroSim
    print ("\n----Calling NeuroSim----\n")
    drift = ['./neurosim.sh']

    if var == 0:
        drift.append(drift_file(0,0).strip())
    else:
        for i in range(len(vals)):
            drift.append(drift_file(var,vals[i]).strip())

    filename = "AccuracyOutput/console/" + str(var) + ".txt"
    file_=open(filename, 'w+')
    subprocess.check_call(drift,stderr=file_)
    fclose(filename)
    return
    
# Main
def main():
    
    arr = np.genfromtxt('run_summary.txt', dtype=np.int32)
    if not arr:
        callNeuroSim(0,0)
    else:
        callNeuroSim(arr[0],arr[1:])

    return
        
if __name__ == "__main__":
    main()