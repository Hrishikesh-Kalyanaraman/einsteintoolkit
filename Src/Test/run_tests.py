#!/usr/bin/env python3

import sys
import os
import re
import shutil

# This function traverses the "Test" directory and returns a list of 
# directory names. One for each test to run.
def find_tests():
    #Check if the "Test" directory exists.
    if not os.path.isdir("Test"):
        print("Error: the directory 'Test' does not exist")
        print("The script has to run from the SelfForce-1D root directory")
        sys.exit(-1)
    # Initialize empty list.
    dirs = []
    # Join the root directory name and the directory name to form the list
    # of directories.
    for (dirpath, dirnames, filenames) in os.walk("Test"):
        for name in dirnames:
            dirs.append(os.path.join(dirpath,name))

    return dirs


# This function checks for the 'pattern' in the 'filename' and adds it to the
# 'dictionary' with 'dirname' as the key if it is found. It throws an error
# and exits if such a match has been found previously for this 'dirname'
def match_pattern(dictionary, pattern, dirname, filename, err_string):
    # Search for the pattern in the filename.
    res = re.search(pattern, filename)
    # If there is a match...
    if ( res != None ):
        # Check if that directory already has an entry in the dictionary.
        if ( dirname in dictionary ):
            # If there is already an entry, throw an error and exit as there
            # should only be one file of this type.
            print("Error: multiple "+err_string+" files in "+dirname)
            sys.exit(1)
        else:
            # If there isn't an entry, create one.
            dictionary[dirname] = os.path.join(dirname,filename)


# This function checks for the 'pattern' in the 'filename' and adds it to the
# 'dictionary' with 'dirname' as the key if it is found. The 'filename' are
# stored as a list in the dictionary so multiple files matching the 'pattern'
# are allowed.
def find_data_files(dictionary, pattern, dirname, filename):
    # Search for the pattern in the filename.
    res = re.search(pattern, filename)
    # If there is a match...
    if ( res != None ):
        # Check if that directory already has an entry in the dictionary.
        if ( dirname in dictionary ):
            # If there is an entry, append the filename to the list.
            dictionary[dirname].append(filename)
        else:
            # If there isn't an entry, create an one element list.
            dictionary[dirname] = [filename]


# This function creates dictionaries for the different kinds of files 
# that can exist in a test directory. It locates the parameter file (.par),
# the test definition file (.test) and the data files (.asc), creates
# dictionaries for each filetype and returns the dictionaries as a 3-tuple.
# Multiple parameter or test definition files are not allowed as the testing
# would not know which one to use. Other file types are currently ignored.
def find_files(dirs):
    # Create empty dictionaries.
    par = {}
    test = {}
    data = {}
    inputs = {}
    # The patterns for the files to search for.
    pattern1 = r"\w*\.par"
    pattern2 = r"\w*\.test"
    pattern3 = r"\w*\.asc"
    pattern4 = r"\w*\.dat"
    # Loop over the directories
    for dirname in dirs:
        # Search for parameter, test description and data files and store them
        # in dictionaries with the directory name as the key.
        for (dirpath, dirnames, filenames) in os.walk(dirname):
            for filename in filenames: 
                # Does the name match the parameter file?
                match_pattern ( par, pattern1, dirname, filename, "parameter" )
                # Does the name match the test definition file?
                match_pattern ( test, pattern2, dirname, filename, "test definition" )
                # Does the name match a data file?
                find_data_files ( data, pattern3, dirname, filename )
                # Does the name match an input file?
                find_data_files ( inputs, pattern4, dirname, filename )
        # If 'dirname' is not in the inputs dictionary...
        if not dirname in inputs:
            # Create an entry with an empty list.
            inputs[dirname] = []
        # If 'dirname' is not in the data dictionary...
        if not dirname in data:
            # Print an error message.
            print("Error: "+dirname+" has no data files")
            # And abort.
            sys.exit(-1)

    return (par,test,data,inputs)


# This function reads the test definition file and extracts the information
# about the executable name and the absolute and relative tolerences to use
# for this test and returns a 3-tuple with that information.
# The function needs to have some checks and error handling added.
def parse_test_definitions(dirs,testfiles):
    # Create empty dictionaries.`
    exe = {}
    relerr = {}
    abserr = {}
    # Loop over the directories.
    for dirname in dirs:
#        print("Opening "+testfiles[dirname])
        # Open the test definition file if it exists, otherwise abort.
        if dirname in testfiles:
            file_object = open(testfiles[dirname],"r")
        else:
            print("Test definition file missing for "+dirname)
            sys.exit(-1)
        # Read all the lines in the file and strip of the newline from the end
        # of each line.
        lines = [line.rstrip() for line in file_object.readlines()]
        # Close the file.
        file_object.close()
        # Loop over the lines.
        for line in lines:
            # Split the line into words using ":\t" as the separator.
            words = line.split(":\t")
#            print(words)
            # Find the line that defines the executable to use.
            if ( words[0]=="Executable" ):
                exe[dirname] = words[1]
            # Find the line that defines the relative tolerance to use.
            if ( words[0]=="Relative tolerance" ):
                relerr[dirname] = float(words[1])
            # Find the line that defines the absolute tolerance to use.
            if ( words[0]=="Absolute tolerance" ):
                abserr[dirname] = float(words[1])
    return ( exe, relerr, abserr ) 


# This function runs the test in 'dirname' using 'exe' as the executable,
# 'par' as the parameter file and then compares the resulting files with
# the data files in 'data' using 'relerr' and 'abserr' as relative and
# absolute tolerances. It reports on whether the test passes or not with
# some information about where the failures are located.
def run_test(dirname,exe,par,data,inputs,relerr,abserr):
    separator = "/"
    # Split the full path to the parameter file.
    par_split = par.split(separator)
    # Copy the split parameter file to the split destination.
    destination_split=par_split
    # Replace the root directory with "RunTests" as we want to run the test
    # there.
    destination_split[0] = "RunTests"
    # Join back the path to the destination directory.
    destination = separator.join(destination_split[:-1])
    print(destination)
    # If the destination directory exists, remove it.
    if os.path.exists(destination):
        shutil.rmtree(destination)
    # Create the destination directory.
    os.makedirs(destination,exist_ok=True)
    # Copy the parameter file to the destination directory.
    shutil.copy2(par,destination)
    # Copy any input files.
    copy_input_files(inputs,dirname,destination)
    # Create the relative path (from the destination directory) to the
    # executable.
    executable = "../../Build/Exe/"+exe
    # Change to the destination directory.
    os.chdir(destination)
    # Create the command to run the executable.
    command = executable+" "+par_split[-1]+" >& "+re.sub(r'par',
                                                         r'out',par_split[-1])
    print("Running "+dirname+":")
    # As we are redirecting stdout and stderr to file anyway we just use
    # os.system to run the command.
    os.system("bash -c '%s'" % command)
    relpath_dir = "../../"+dirname
    # Compare the test data.
    (fail,tol,failures) = compare_data(relpath_dir,data,relerr,abserr)
    # Check for any failures but distinguish between identical results and
    # passing within tolerances. A missing file counts as a failure.
    # Check if all files pass.
    if fail==0:
        # Check if all tests passes with identical results.
        if tol==0:
            print("Passed: "+str(len(data))+" files identical")
        # Otherwise report on how many files pass only to within tolerances.
        else:
            print("Passed: differences below tolerance in "+str(tol)+" files")
    else:
        print("Failed: differences above tolerance in "+str(fail)+" files")
    # If there is any failure information to report...
    if (len(failures)>0):
        # Loop over the failure messages.
        for failure in failures:
            # Print the failure messages.
            print(failure)
    # Print an empty line between tests.
    print("")
    # Change back to the root directory.
    os.chdir("../..")


# This function copies all the input data files in 'inputs' from 'dirname' to
# 'destination'.
def copy_input_files(inputs,dirname,destination):
    # Loop over the filenames
    for filename in inputs:
        # Join the data input file name with the directory name.
        source = os.path.join(dirname,filename)
        # Copy the file.
        shutil.copy2(source,destination)


# This function loops over all filenames listed in 'data' and compares
# them with the corresponding file in 'path' (it's a relative path) using
# 'relerr' and 'abserr' as the relative and absolute tolerance.
def compare_data(path,data,relerr,abserr):
    # Create and empty list of failures.
    failures = []
    # Initialize a counter for failed files to zero.
    filefail = 0
    # Initialize a counter for files that pass only to within thentolerances
    # to zero.
    filetol = 0
    # Loop over the files present in the test directory.
    for filename in data:
        # Check that the corresponding file exists in the current directory.
        if os.path.exists(filename) and os.path.isfile(filename):
            # Compare the files and receive information about fail/pass/tol.
            (fail,tol) = compare_files(os.path.join(path,filename),
                                       filename,relerr,abserr,failures)
            # If there is a failure increment filefail.
            if fail>0:
                filefail += 1
            # If the file only passes to within tolerances increment filetol.
            if tol>0:
                filetol += 1
        else:
            # Append an error message to failures with the information about
            # this missing file
            failures.append("Data file "+filename+" is missing")
            # Increment filefail.
            filefail += 1

    # Return information in a 3-tuple.
    return (filefail,filetol,failures)


# This function campares the numbers within 2 files ('base' and 'new') to
# within relative tolerance 'relerr' and absolute tolerance 'abserr' and
# updates the list of 'failures'
def compare_files(base,new,relerr,abserr,failures):
    # Open the reference file.
    basefile = open(base,"r")
    # Open the newly created file.
    newfile = open(new,"r")
    # Initialize counters to zero for linenumber (the line number), failcount
    # (the number of lines with failures) and tolcount (the number of lines
    # that pass only to tolerences).
    linenumber = 0
    failcount = 0
    tolcount = 0
    # A regex to test if the first item on a line is a number (it ignores
    # exponents).
    numberregex = re.compile("^\s*-?\d+\.?\d*")
    # Loop over all lines in the basefile.
    for line in basefile:
        # increment the line number.
        linenumber += 1
        # Read the corresponding lines in the new file.
        newline = newfile.readline()
        # If there is a line comare it with the line from the base file.
        if newline:
            # Search for a number in the beginning of the line.
            res = numberregex.search(line)            
            # If there is a number, compare the line number by number
            if ( res != None ):
                # Create a list of the numbers in the line in the base file.
                val_base = [float(x) for x in line.split()]
                # Create a list of the numbers in the line in the new file.
                val_new = [float(x) for x in newline.split()]
                # Compare the list of numbers.
                (fail,tol) = compare_numbers(val_base,val_new,new,linenumber,
                                             relerr,abserr,failures)
                # If there is a failure...
                if fail:
                    # increment failcount.
                    failcount += 1
                # If a number only passes to tolerance...
                if tol:
                    # increment tolcount.
                    tolcount += 1
        else:
            # If the new file has less lines than the base file add a failure
            # message.
            failures.append(new+" has less lines than "+base)
            # And break. We are done with this file.
            break

    # If there is one failure.
    if failcount==1:
        # Append a failure message.
        failures.append(new+" has "+str(failcount)+" difference above tolerance")
    # If there are multiple failures.
    if failcount>1:
        # Append a failure message.
        failures.append(new+" has "+str(failcount)+" differences above tolerance")
    # If there is one line that passes only to within tolerances.
    if tolcount==1:
        # Append a nonfailure message.
        failures.append(new+" has "+str(tolcount)+" difference belov tolarance")
    # If there are multiple lines that passes only to within tolerances.
    if tolcount>1:
        # Append a nonfailure message.
        failures.append(new+" has "+str(tolcount)+" differences belov tolarance")

    # return the counts of failures and tolerance passes.
    return (failcount, tolcount)


# This function compares the elements of two lists of numbers ('x' and 'y')
# coming from line number 'n' in 'filename' within relative tolerance 'reltol'
# and absolute tolerance 'abstol'. Failures get appended to 'failures'.
def compare_numbers(x,y,filename,n,reltol,abstol,failures):
    # Initialize the boolean variables 'fail' and 'tol'
    fail = False
    tol = False
    # Get the length of 'x'.
    lenx = len(x)
    # Get the length of 'y'.
    leny = len(y)
    # If there are diffent number of elements in 'x' and 'y'
    if lenx!=leny:
        # Append a failure message.
        failures.append(filename+": different number of values at line "+str(n))
        # And indicate a failure.
        fail = True
    else:
        # Compare the numbers element by element.
        for i in range(lenx):
            # Calculate the absolute difference.
            abserr = abs(y[i]-x[i])
            # Initialize relerr to zero in case it's not calculated later.
            relerr = 0.0
            # If the absolute difference is above tolerance...
            if abserr>abstol:
                # As the absolute difference could come from comparing 2 large
                # numbers we then also calculate the relative difference using
                # division with a number that is not zero.
                if x[i]>0:
                    relerr = abs(y[i]/x[i]-1.0)
                else:
                    relerr = abs(x[i]/y[i]-1.0)
               
                # If the relative difference is also larger than the
                # tolerance...
                if relerr>reltol:
                   # Set fail to True.
                   fail = True
            # If the absolute error is not zero and either the absolute 
            # difference or relative difference are below tolerances...       
            if (abserr>0.0 and (abserr<abstol or relerr<reltol)):
                # The comparison passes to within tolerances...
                tol = True
      
    # If there is a failure...
    if fail:
        # Append a failure message indicating where the failure occurred.
        failures.append(filename+": values differ significantly at line "+str(n))
    # Return the information about the comparison.
    return (fail,tol)


# Find the tests to run.
dirs = find_tests()

# Find the parameter, test definition and data files for all the tests.
(parfiles,testfiles,datafiles,inputfiles) = find_files(dirs)

# Parse the test definition file.
(executables,relative_error,absolute_error) = parse_test_definitions(
    dirs,testfiles)

# Print information about the number of tests.
print("Running "+str(len(dirs))+" tests\n")

# Loop over the tests.
for dirname in dirs:
    # Run the tests.
    run_test(dirname,executables[dirname],parfiles[dirname],datafiles[dirname],
             inputfiles[dirname],relative_error[dirname],
             absolute_error[dirname])
