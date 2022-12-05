# Import modules
import sys
import csv
import glob
import matlab.engine
import numpy as np

# Update csv file value
def writefile(var, val):
    """
    var - what variable user is sweeping
    val - single value to change in csv files
    
    Goes through all chip files in inputs/chips to write over csv files
    """
    # Go through all files
    files = glob.glob("TSV_stacked_3D_v4/air_cooling/inputs/chips/*.csv")
    for file in files:
        with open(file, 'rt') as f:
            rows = csv.reader(f, skipinitialspace=False, delimiter=',', quoting=csv.QUOTE_NONE)
            data = [data for data in rows]
            
            if var == 1:
                data[7][3] = val * 10E-6 # chip thickness
            elif var == 2:
                data[13][3] = val * 10E-3 # X-size
                data[14][3] = val * 10E-3 # Y-size
            elif var == 3:
                data[22][3] = val * 10E-6 # TSV Diameter
            elif var == 4:
                # TSV Pitch
                data[24][3] = val * 10E-6 # X
                data[25][3] = val * 10E-6 # Y
            elif var == 5:
                data[16][3] = val * 10E-6 # Bump Diameter
            else:
                # Bump Pitch
                data[17][3] = val * 10E-6 # X
                data[18][3] = val * 10E-6 # Y
        
        # Write value to each file
        temp = csv.writer(open(file, 'wt', newline=''))
        temp.writerows(data)
       
    return

def matlabengine(var, vals):
    """
    var - what variable user is sweeping
    vals - array of values to sweep
    
    Calls matlab engine to get max junction temperature and retention data
    """
    
    # Call MATLAB engine and edit file values
    # Single run
    if var == 0:
        print ("***CALLING MATLAB ENGINE***")
        print ("------Single Run------")
        eng = matlab.engine.start_matlab()
        eng.startup(var, 0, nargout=0)
        print ("***MATLAB ENGINE HAS EXITED***")
        eng.quit()
        
    else:
        print ("***CALLING MATLAB ENGINE***")
        for i in range(len(vals)):
            # Write to csv file
            writefile(var, vals[i])
            
            # Chip thickness
            if var == 1:
                print ("----Chip Thickness: " + str(vals[i]) + "um----\n")
                # Call matlab engine
                eng = matlab.engine.start_matlab()
                eng.startup(var, vals[i], nargout=0)
                eng.quit()
            elif var == 2:
                # Chip Size
                print ("----Chip Size: " + str(vals[i]) + "mm----\n")
                # Call matlab engine
                eng = matlab.engine.start_matlab()
                eng.startup(var, vals[i], nargout=0)
                eng.quit()
            elif var == 3:
                # TSV Diameter
                print ("----TSV Diameter: " + str(vals[i]) + "um----\n")
                # Call matlab engine
                eng = matlab.engine.start_matlab()
                eng.startup(var, vals[i], nargout=0)
                eng.quit()
            elif var == 4:
                # TSV Pitch
                print ("----TSV Pitch: " + str(vals[i]) + "um----\n")
                # Call matlab engine
                eng = matlab.engine.start_matlab()
                eng.startup(var, vals[i], nargout=0)
                eng.quit()
            elif var == 5:
                # Bump Diameter
                print ("----Bump Diameter: " + str(vals[i]) + "um----\n")
                # Call matlab engine
                eng = matlab.engine.start_matlab()
                eng.startup(var, vals[i], nargout=0)
                eng.quit()
            else:
                # Bump Pitch
                print ("----Bump Pitch: " + str(vals[i]) + "um----\n")
                # Call matlab engine
                eng = matlab.engine.start_matlab()
                eng.startup(var, vals[i], nargout=0)
                eng.quit()

        print ("***MATLAB ENGINE HAS EXITED***")
        
    return

# Main
def main():
    
    # Usage statement
    if len(sys.argv) < 2:
        print ("Usage: NAME VAR val1 val2 ...")
        print ("\nVAR: User design space exploration")
        print ("\n\t 0 - single run \n\t 1 - chip thickness: X(um) \n\t 2 - chip size: XY(mm)")
        print ("\t 3 - TSV Diameter: X(um) \n\t 4 - TSV Pitch: X(um) \n\t 5 - Bump Diameter: X(um) \n\t 6 - Bump Pitch: X(um)")
        return

    var = int(sys.argv[1])

    # Check for invalid entries
    if (var > 6 or var < 0):
        print ("ERROR: Invalid entry.")
        return
    
    # Get range of values, check for invalid entry
    vals = [float(i) for i in sys.argv[2:]]
    if (not vals and var):
        print ("ERROR: Missing values. Please add arguments.")
        return

    # Calling matlab engine
    matlabengine(var, vals)
    np.savetxt('run_summary.txt', [int(i) for i in sys.argv[1:]], fmt='%1.0f')

    return
        
if __name__ == "__main__":
    main()