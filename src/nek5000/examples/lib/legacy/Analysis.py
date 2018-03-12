#! /usr/bin/python
# Python module to run top-down tests for Nek

import sys
import os

###############################################################################
def Test(name, logfile,listOfValue)  :
    """A Test function which look in the log file and compare the value to the target value 

        --Variable :
            name (string): name of the test
            logfile (string) : path of the log file
            listOfValue (list) : list of the different ['name',target,tolerance] we want to check
         --Function :
           Test will read the logfile and for all the set of ['name',target,tolerance,position] in listOfValue
           will catch the value of 'name', compare it to target according to tolerance and return
           success or failure.
           Position is the number of the row (starting from the right in which the information is."""
           
    global num_test
    global num_success
    
    test_result = False
    reported_ValueError = False
    reported_IndexError = False
    numTest = len(listOfValue)                             #Number of tests to do
    success = 0
    num_test += numTest
    #Test that the input file are here       
    try :
        olog = open(logfile, 'r')
        olog.close()
        inputisgood = True
    except IOError :
        print("[%s]...Sorry, I must skip this test."%name)
        print("[%s]...The logfile is missing or doesn't have the correct name..."%name)
        inputisgood = False

    #If we could find the log file, then we run the test    
    if inputisgood :        
        log = open(logfile,'r')
        for line in log :                                               #we read the log file line by line  
            for set in listOfValue :                                    #for each line we loof if the value we want are there
                if set[0] in line :                                     #set[0] is the name of the value
                    # Try to parse testvalue.  If there is a parsing error in
                    # line, assume that set wasn't found in line and continue.
                    try:
                        testvalue = float(line.split()[-set[3]])                 #the value is on the set[3]th column from the right of the line 
                    except (ValueError):
                        if not reported_ValueError:
                            print("Warning: Attempted to parse non-numerical value for test \"%s\".  Logfile may be malformatted"%name)
                            reported_ValueError = True
                    except (IndexError):
                        if not reported_IndexError:
                            print("Warning: Fewer columns than excpected for test \"%s\".  Logfile may be malformatted"%name)
                            reported_IndexError = True
                    else:
                        print("[%s] %s : %s"%(name,set[0],testvalue))
                        if (abs(testvalue - set[1]) < set[2]) :             #set[1] is the target value / set[2] is the tolerance
                           if (testvalue != 0.0) :                          #Checks that it is not 0.0(failure)
                              success += 1
                              num_success +=1
                        listOfValue.remove(set)                             #The value was found so we remove it from the search
        log.close()
    if success == numTest :
        test_result = True
    elif (len(listOfValue) > 0) :
        print("[%s]...I couldn't find all the requested value in the log file..."%name)
        slist = ""
        for set in listOfValue :
            slist = slist + set[0] + ", "
        print("[%s]...%s were not found..."%(name,slist))

    return test_result
    
###############################################################################
def Run(name, logfile,listOfValue)  :
    """A Run function which runs the test and prints out the result
        --Variable :
            name (string): name of the test
            logfile (string) : path of the log file
            listOfValue (list) : list of the different ['name',target,tolerance] we want to check"""
           
    Result = Test(name,logfile, listOfValue)   
    if Result :
        print("%s : ."%name)
    else :
        print("%s : F "%name)
###############################################################################
def FindPhrase(name, logfile, keyword) :
    """A  Test to search the logfile for a specific word or phrase
        --Variable :
            name (string): name of the test
            logfile (string): path of the logfile
            keyword (string): word or phrase searching for
        --Function :
            Tests to find a phrase or word in the logfile"""

    global num_test
    global num_success

    num_test += 1
    result    = False

    #Test that the input file is here
    try :
        olog = open(logfile, 'r')
        olog.close()
        inputisgood = True
    except IOError :
        print("[%s]...Sorry, I must skip this test."%name)
        print("[%s]...The logfile is missing or doesn't have the correct name..."%name)
        inputisgood = False

    #If we could find the log file, then we run the test
    if inputisgood :
        openlog  = open(logfile,'r')
        for line in openlog :
            if keyword in line :
                num_success += 1
                result       = True
                print("[%s] : %s"%(name,keyword))
        openlog.close()

        if result :
            print("%s : ."%name)                 #prints the result
        else :
            print("[%s]...I couldn't find '%s' in the logfile..."%(name,keyword))
            print("%s : F "%name)

###############################################################################
def DFdPhrase(name, logfile, keyword) :
    """A  Test to search the logfile for a specific word or phrase, returns True, if not found
        --Variable :
            name (string): name of the test
            logfile (string): path of the logfile
            keyword (string): word or phrase searching for
        --Function :
            Tests to find that a phrase or word is not the logfile"""

    global num_test
    global num_success

    num_test += 1
    result    = False

    #Test that the input file is here
    try :
        olog = open(logfile, 'r')
        olog.close()
        inputisgood = True
    except IOError :
        print("[%s]...Sorry, I must skip this test."%name)
        print("[%s]...The logfile is missing or doesn't have the correct name..."%name)
        inputisgood = False

    #If we could find the log file, then we run the test
    if inputisgood :
        openlog  = open(logfile,'r')
        for line in openlog :
            if keyword in line :
                result       = True
                print("[%s] : %s"%(name,keyword))
        openlog.close()

        if result :
            print("%s : F"%name)                 #prints the result
        else :
            num_success += 1
            print("%s : . "%name)

###############################################################################
###############################################################################
#  Check if mpi tests were run..
TestsToDo = sys.argv
if "mpi" in TestsToDo:
   ifmpi = True 
else:
   ifmpi = False
   print("NO MPI TESTS BEING RAN! ")
   print("If incorrect, call Analysis with 'mpi' as an argument")
###############################################################################
num_test = 0 
num_success = 0
print("Beginning of top-down testing\n\n")
print("    . : successful test, F : failed test\n\n")
###############################################################################
#---------------------Tools----------------------------------------------------
#---------------------Tests----------------------------------------------------
print("\nBEGIN TESTING TOOLS")  
#SRL Compiler
log = "./tools.out"
value = "Error "
DFdPhrase("Tools",log,value)

###############################################################################
#---------------------MPI------------------------------------------------------
#---------------------Pn-Pn----------------------------------------------------
print("\n\n2d_eig Example")  
#SRL
log = "./srlLog/eig1.err"
value = [[' 2   ',0,21,6],
         [' 3   ',0,39,6],
         [' 6  ' ,0,128,6],
         [' 10  ',0,402,6]]
Run("Example 2d_eig/SRL: Serial-iter/err",log,value)
value = [[' 2   ',1.6329577E-01,1e-01,2],
         [' 3   ',1.4437432E-01,1e-01,2],
         [' 6  ' ,1.1637760E-02,1e-02,2],
         [' 10  ',6.8264260E-06,1e-06,2]]
Run("Example 2d_eig/SRL: Serial-iter/err",log,value)

#SRL2
log = "./srl2Log/eig1.err"
value = [[' 2   ',0,21,6],
         [' 3   ',0,39,6],
         [' 6  ' ,0,126,6],
         [' 10  ',0,402,6]]
Run("Example 2d_eig/SRL2: Serial-iter/err",log,value)

value = [[' 2   ',1.6329577E-01,1e-01,2],
         [' 3   ',1.4437432E-01,1e-01,2],
         [' 6  ' ,1.1637760E-02,1e-02,2],
         [' 10  ',6.8264260E-06,1e-06,2]]
Run("Example 2d_eig/SRL2: Serial-iter/err",log,value)




print("\n\n3Dbox Example")  
#MPI
if ifmpi:

    log = "./mpiLog/b3d.log.1"
    value = "end of time-step loop"
    FindPhrase("Example 3dbox/MPI: Serial",log,value)

    log = "./mpiLog/b3d.log.4"
    value = "end of time-step loop"
    FindPhrase("Example 3dbox/MPI: Parallel",log,value)

#SRL
log = "./srlLog/b3d.log.1"
value = "end of time-step loop"
FindPhrase("Example 3dbox/SRL: Serial",log,value)

#MPI2
if ifmpi:
    log = "./mpi2Log/b3d.log.1"
    value = "end of time-step loop"
    FindPhrase("Example 3dbox/MPI2: Serial",log,value)

    log = "./mpi2Log/b3d.log.4"
    value = "end of time-step loop"
    FindPhrase("Example 3dbox/MPI2: Parallel",log,value)

#SRL2
log = "./srl2Log/b3d.log.1"
value = "end of time-step loop"
FindPhrase("Example 3dbox/SRL2: Serial",log,value)




print("\n\naxi Example")  
#MPI
if ifmpi:

    log = "./mpiLog/axi.log.1"
    value = [['PRES: ',0,76,4]]
    Run("Example axi/MPI: Serial-iter",log,value)

    log = "./mpiLog/axi.log.4"
    value = [['PRES: ',0,76,4]]
    Run("Example axi/MPI: Parallel-iter",log,value)

#SRL
log = "./srlLog/axi.log.1"
value = [['total solver time',0.1,2,2],
         ['PRES: ',0,76,4]]
Run("Example axi/SRL: Serial-time/iter",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/axi.log.1"
    value = [['U-Press ',0,104,5]]
    Run("Example axi/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/axi.log.4"
    value = [['U-Press ',0,104,5]]
    Run("Example axi/MPI2: Parallel-iter",log,value)

#SRL2
log = "./srl2Log/axi.log.1"
value = [['total solver time',0.1,4,2],
         ['U-Press ',0,104,5]]
Run("Example axi/SRL2: Serial-iter",log,value)




print("\n\nbenard-ray_9 Example")  
#MPI
if ifmpi:

    log = "./mpiLog/ray_9.log.1"
    value = [['gmres: ',0,23,7]]
    Run("Example benard/ray_9/MPI: Serial-iter",log,value)

#SRL
log = "./srlLog/ray_9.log.1"
value = [['total solver time',0.1,30,2],
         ['gmres: ',0,23,7]]
Run("Example benard/ray_9/SRL: Serial-iter",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/ray_9.log.1"
    value = [['gmres: ',0,11,6]]
    Run("Example benard/ray_9/MPI2: Serial-iter",log,value)

#SRL2
log = "./srl2Log/ray_9.log.1"
value = [['total solver time',0.1,40,2],
         ['gmres: ',0,11,6]]
Run("Example benard/ray_9/SRL2: Serial-time/iter",log,value)




print("\n\nbenard-ray_dd Example")  
#MPI
if ifmpi:

    log = "./mpiLog/ray_dd.log.1"
    value = [['gmres: ',0,11,6]]
    Run("Example benard/ray_dd/MPI: Serial-iter",log,value)

    log = "./mpiLog/benard.err"
    value = [['ray_dd.log.1',1707.760,1,7]]
    Run("Example benard/ray_dd/MPI: Serial-error",log,value)

#SRL
log = "./srlLog/ray_dd.log.1"
value = [['total solver time',0.1,24,2],
         ['gmres: ',0,11,6]]
Run("Example benard/ray_dd/SRL: Serial-time/iter",log,value)

log = "./srlLog/benard.err"
value = [['ray_dd.log.1',1707.760,1,7]]
Run("Example benard/ray_dd/SRL: Serial-error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/ray_dd.log.1"
    value = [['gmres: ',0,11,6]]
    Run("Example benard/ray_dd/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/benard.err"
    value = [['ray_dd.log.1',1707.760,1,7]]
    Run("Example benard/ray_dd/MPI2: Serial-error",log,value)

#SRL2
log = "./srl2Log/ray_dd.log.1"
value = [['total solver time',0.1,20,2],
         ['gmres: ',0,11,6]]
Run("Example benard/ray_dd/SRL2: Serial-time/iter",log,value)

log = "./srl2Log/benard.err"
value = [['ray_dd.log.1',1707.760,1,7]]
Run("Example benard/ray_dd/SRL2: Serial-error",log,value)




print("\n\nbenard-ray_dn Example")  
#MPI
if ifmpi:

    log = "./mpiLog/ray_dn.log.1"
    value = [['gmres: ',0,11,6]]
    Run("Example benard/ray_dn/MPI: Serial-iter",log,value)

    log = "./mpiLog/benard.err"
    value = [['ray_dn.log.1',1100.650,1,7]]
    Run("Example benard/ray_dn/MPI: Serial-error",log,value)

#SRL
log = "./srlLog/ray_dn.log.1"
value = [['total solver time',0.1,30,2],
         ['gmres: ',0,11,6]]
Run("Example benard/ray_dn/SRL: Serial-time/iter",log,value)

log = "./srlLog/benard.err"
value = [['ray_dn.log.1',1100.650,1,7]]
Run("Example benard/ray_dn/SRL: Serial-error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/ray_dn.log.1"
    value = [['gmres: ',0,11,6]]
    Run("Example benard/ray_dn/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/benard.err"
    value = [['ray_dn.log.1',1100.650,1,7]]
    Run("Example benard/ray_dn/MPI2: Serial-error",log,value)

#SRL2
log = "./srl2Log/ray_dn.log.1"
value = [['total solver time',0.1,12,2],
         ['gmres: ',0,11,6]]
Run("Example benard/ray_dn/SRL2: Serial-time/iter",log,value)

log = "./srl2Log/benard.err"
value = [['ray_dn.log.1',1100.650,1,7]]
Run("Example benard/ray_dn/SRL2: Serial-error",log,value)




print("\n\nbenard-ray_nn Example")  
#MPI
if ifmpi:

    log = "./mpiLog/ray_nn.log.1"
    value = [['gmres: ',0,14,6]]
    Run("Example benard/ray_nn/MPI: Serial-iter",log,value)

    log = "./mpiLog/benard.err"
    value = [['ray_nn.log.1',657.511,.1,7]]
    Run("Example benard/ray_nn/MPI: Serial-error",log,value)

#SRL
log = "./srlLog/ray_nn.log.1"
value = [['total solver time',0.1,30,2],
         ['gmres: ',0,14,6]]
Run("Example benard/ray_nn/SRL: Serial-time/iter",log,value)

log = "./srlLog/benard.err"
value = [['ray_nn.log.1',657.511,.1,7]]
Run("Example benard/ray_nn/SRL: Serial-error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/ray_nn.log.1"
    value = [['gmres: ',0,14,6]]
    Run("Example benard/ray_nn/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/benard.err"
    value = [['ray_nn.log.1',657.511,.1,7]]
    Run("Example benard/ray_nn/MPI2: Serial-error",log,value)

#SRL2
log = "./srl2Log/ray_nn.log.1"
value = [['total solver time',0.1,20,2],
         ['gmres: ',0,14,6]]
Run("Example benard/ray_nn/SRL2: Serial-time/iter",log,value)

log = "./srl2Log/benard.err"
value = [['ray_nn.log.1',657.511,.1,7]]
Run("Example benard/ray_nn/SRL2: Serial-error",log,value)




print("\n\nblasius Example")  
#MPI
if ifmpi:

    log = "./mpiLog/blasius.log.1"
    value = [['gmres: ',0,162,7]]
    Run("Example blasius/MPI: Serial-iter",log,value)

    log = "./mpiLog/blasius.err.1"
    value = [['delta',1.26104E+00,1e-05,5]]
    Run("Example blasius/MPI: Serial-error",log,value)

    log = "./mpiLog/blasius.log.4"
    value = [['gmres: ',0,162,7]]
    Run("Example blasius/MPI: Parallel-iter",log,value)

    log = "./mpiLog/blasius.err.4"
    value = [['delta',1.26104E+00,1e-05,5]]
    Run("Example blasius/MPI: Parallel-error",log,value)

#SRL
log = "./srlLog/blasius.log.1"
value = [['total solver time',0.1,30,2],
         ['gmres: ',0,162,7]]
Run("Example blasius/SRL: Serial-time/iter",log,value)

log = "./srlLog/blasius.err.1"
value = [['delta',1.26104E+00,1e-05,5]]
Run("Example blasius/SRL: Serial-error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/blasius.log.1"
    value = [['gmres: ',0,125,6]]
    Run("Example blasius/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/blasius.err.1"
    value = [['delta',1.26104E+00,1e-05,5]]
    Run("Example blasius/MPI2: Serial-error",log,value)

    log = "./mpi2Log/blasius.log.4"
    value = [['gmres: ',0,125,6]]
    Run("Example blasius/MPI2: Parallel-iter",log,value)

    log = "./mpi2Log/blasius.err.4"
    value = [['delta',1.26104E+00,1e-05,5]]
    Run("Example blasius/MPI2: Parallel-error",log,value)

#SRL
log = "./srl2Log/blasius.log.1"
value = [['total solver time',0.1,30,2],
         ['gmres: ',0,125,6]]
Run("Example blasius/SRL2: Serial-time/iter",log,value)

log = "./srl2Log/blasius.err.1"
value = [['delta',1.26104E+00,1e-05,5]]
Run("Example blasius/SRL2: Serial-error",log,value)




print("\n\nconj_ht Example")  
#MPI
if ifmpi:

    log = "./mpiLog/conj_ht.log.1"
    value = [['gmres: ',0,46,7]]
    Run("Example conj_ht/MPI: Serial-iter",log,value)

    log = "./mpiLog/conj_ht.err.1"
    value = [['tmax',1.31190E+01,1e-06,2]]
    Run("Example conj_ht/MPI: Serial-error",log,value)

    log = "./mpiLog/conj_ht.log.4"
    value = [['gmres: ',0,46,7]]
    Run("Example conj_ht/MPI: Parallel-iter",log,value)

    log = "./mpiLog/conj_ht.err.4"
    value = [['tmax',1.31190E+01,1e-06,2]]
    Run("Example conj_ht/MPI: Parallel-error",log,value)

#SRL
log = "./srlLog/conj_ht.log.1"
value = [['total solver time',0.1,7,2],
         ['gmres: ',0,46,7]]
Run("Example conj_ht/SRL: Serial-time/iter",log,value)

log = "./srlLog/conj_ht.err.1"
value = [['tmax',1.31190E+01,1e-06,2]]
Run("Example conj_ht/SRL: Serial-error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/conj_ht.log.1"
    value = [['gmres: ',0,26,6]]
    Run("Example conj_ht/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/conj_ht.err.1"
    value = [['tmax',1.31190E+01,1e-06,2]]
    Run("Example conj_ht/MPI2: Serial-error",log,value)

    log = "./mpi2Log/conj_ht.log.4"
    value = [['gmres: ',0,26,6]]
    Run("Example conj_ht/MPI2: Parallel-iter",log,value)

    log = "./mpi2Log/conj_ht.err.4"
    value = [['tmax',1.31190E+01,1e-06,2]]
    Run("Example conj_ht/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/conj_ht.log.1"
value = [['total solver time',0.1,7,2],
         ['gmres: ',0,26,6]]
Run("Example conj_ht/SRL2: Serial-time/iter",log,value)

log = "./srl2Log/conj_ht.err.1"
value = [['tmax',1.31190E+01,1e-06,2]]
Run("Example conj_ht/SRL2: Serial-error",log,value)




print("\n\ncone016 Example")  
#MPI
if ifmpi:

    log = "./mpiLog/cone016.err.1"
    value = [['Tmax',8.5065E-01,1e-06,2]]
    Run("Example cone016/MPI: Serial-error",log,value)

    log = "./mpiLog/cone016.err.4"
    value = [['Tmax',8.5065E-01,1e-06,2]]
    Run("Example cone016/MPI: Parallel-error",log,value)

#SRL
log = "./srlLog/cone016.log.1"
value = [['total solver time',0.1,7,2]]
Run("Example cone016/SRL: Serial-time",log,value)

log = "./srlLog/cone016.err.1"
value = [['Tmax',8.5065E-01,1e-06,2]]
Run("Example cone016/SRL: Serial-error",log,value)




print("\n\ncone064 Example")  
#MPI
if ifmpi:

    log = "./mpiLog/cone064.err.1"
    value = [['Tmax',7.9285E-01,1e-06,2]]
    Run("Example cone064/MPI: Serial-error",log,value)

    log = "./mpiLog/cone064.err.4"
    value = [['Tmax',7.9285E-01,1e-06,2]]
    Run("Example cone064/MPI: Parallel-error",log,value)

#SRL
log = "./srlLog/cone064.log.1"
value = [['total solver time',0.1,7,2]]
Run("Example cone064/SRL: Serial-time",log,value)

log = "./srlLog/cone064.err.1"
value = [['Tmax',7.9285E-01,1e-06,2]]
Run("Example cone064/SRL: Serial-error",log,value)




print("\n\ncone256 Example")  
#MPI
if ifmpi:

    log = "./mpiLog/cone256.err.1"
    value = [['Tmax',7.4924E-01,1e-06,2]]
    Run("Example cone256/MPI: Serial-error",log,value)

    log = "./mpiLog/cone256.err.4"
    value = [['Tmax',7.4924E-01,1e-06,2]]
    Run("Example cone256/MPI: Parallel-error",log,value)

#SRL
log = "./srlLog/cone256.log.1"
value = [['total solver time',0.1,9,2]]
Run("Example cone256/SRL: Serial-time",log,value)

log = "./srlLog/cone256.err.1"
value = [['Tmax',7.4924E-01,1e-06,2]]
Run("Example cone256/SRL: Serial-error",log,value)




print("\n\cyl_restart Example")
print("\n\ca:::")
#MPI
if ifmpi:

    log = "./mpiLog/ca.log.1"
    value = [['gmres: ',0,85,7]]
    Run("Example restart-ca/MPI: Serial-iter",log,value)

    log = "./mpiLog/ca.err.1"
    value = [['dragy',5.37986119139E-03,1e-06,4]]
    Run("Example restart-ca/MPI: Serial-error",log,value)

    log = "./mpiLog/ca.log.4"
    value = [['gmres: ',0,85,7]]
    Run("Example restart-ca/MPI: Parallel-iter",log,value)

    log = "./mpiLog/ca.err.4"
    value = [['dragy',5.37986119139E-03,1e-06,4]]
    Run("Example restart-ca/MPI: Parallel-error",log,value)

#SRL
log = "./srlLog/ca.log.1"
value = [['gmres: ',0,85,7]]
Run("Example restart-ca/SRL: Serial-iter",log,value)

log = "./srlLog/ca.err.1"
value = [['dragy',5.37986119139E-03,1e-06,4]]
Run("Example restart-ca/SRL: Serial-error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/ca.log.1"
    value = [['gmres: ',0,29,6]]
    Run("Example restart-ca/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/ca.err.1"
    value = [['dragy',5.09547531705E-02,1e-06,4]]
    Run("Example restart-ca/MPI2: Serial-error",log,value)

    log = "./mpi2Log/ca.log.4"
    value = [['gmres: ',0,29,6]]
    Run("Example restart-ca/MPI2: Parallel-iter",log,value)

    log = "./mpi2Log/ca.err.4"
    value = [['dragy',5.09547531705E-02,1e-06,4]]
    Run("Example restart-ca/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/ca.log.1"
value = [['gmres: ',0,29,6]]
Run("Example restart-ca/SRL2: Serial-iter",log,value)

log = "./srl2Log/ca.err.1"
value = [['dragy',5.09547531705E-02,1e-06,4]]
Run("Example restart-ca/SRL2: Serial-error",log,value)


print("\n\cb:::")
#MPI
if ifmpi:

    log = "./mpiLog/cb.log.1"
    value = [['gmres: ',0,77,7]]
    Run("Example restart-cb/MPI: Serial-iter",log,value)

    log = "./mpiLog/cb.err.1"
    value = [['dragy',5.37986119139E-03,1e-06,4]]
    Run("Example restart-cb/MPI: Serial-error",log,value)

    log = "./mpiLog/cb.log.4"
    value = [['gmres: ',0,77,7]]
    Run("Example restart-cb/MPI: Parallel-iter",log,value)

    log = "./mpiLog/cb.err.4"
    value = [['dragy',5.37986119139E-03,1e-06,4]]
    Run("Example restart-cb/MPI: Parallel-error",log,value)

#SRL
log = "./srlLog/cb.log.1"
value = [['gmres: ',0,77,7]]
Run("Example restart-cb/SRL: Serial-iter",log,value)

log = "./srlLog/cb.err.1"
value = [['dragy',5.37986119139E-03,1e-06,4]]
Run("Example restart-cb/SRL: Serial-error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/cb.log.1"
    value = [['gmres: ',0,28,6]]
    Run("Example restart-cb/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/cb.err.1"
    value = [['dragy',5.09547531705E-02,1e-06,4]]
    Run("Example restart-cb/MPI2: Serial-error",log,value)

    log = "./mpi2Log/cb.log.4"
    value = [['gmres: ',0,28,6]]
    Run("Example restart-cb/MPI2: Parallel-iter",log,value)

    log = "./mpi2Log/cb.err.4"
    value = [['dragy',5.09547531705E-02,1e-06,4]]
    Run("Example restart-cb/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/cb.log.1"
value = [['gmres: ',0,28,6]]
Run("Example restart-cb/SRL2: Serial-iter",log,value)

log = "./srl2Log/cb.err.1"
value = [['dragy',5.09547531705E-02,1e-06,4]]
Run("Example restart-cb/SRL2: Serial-error",log,value)



print("\n\pa:::")
#MPI
if ifmpi:

    log = "./mpiLog/pa.log.1"
    value = [['gmres: ',0,85,7]]
    Run("Example restart-pa/MPI: Serial-iter",log,value)

    log = "./mpiLog/pa.err.1"
    value = [['dragy',5.37986119139E-03,1e-06,4]]
    Run("Example restart-pa/MPI: Serial-error",log,value)

    log = "./mpiLog/pa.log.4"
    value = [['gmres: ',0,85,7]]
    Run("Example restart-pa/MPI: Parallel-iter",log,value)

    log = "./mpiLog/pa.err.4"
    value = [['dragy',5.37986119139E-03,1e-06,4]]
    Run("Example restart-pa/MPI: Parallel-error",log,value)

#SRL
log = "./srlLog/pa.log.1"
value = [['gmres: ',0,85,7]]
Run("Example restart-pa/SRL: Serial-iter",log,value)

log = "./srlLog/pa.err.1"
value = [['dragy',5.37986119139E-03,1e-06,4]]
Run("Example restart-pa/SRL: Serial-error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/pa.log.1"
    value = [['gmres: ',0,29,6]]
    Run("Example restart-pa/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/pa.err.1"
    value = [['dragy',5.09547531705E-02,1e-06,4]]
    Run("Example restart-pa/MPI2: Serial-error",log,value)

    log = "./mpi2Log/pa.log.4"
    value = [['gmres: ',0,29,6]]
    Run("Example restart-pa/MPI2: Parallel-iter",log,value)

    log = "./mpi2Log/pa.err.4"
    value = [['dragy',5.09547531705E-02,1e-06,4]]
    Run("Example restart-pa/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/pa.log.1"
value = [['gmres: ',0,29,6]]
Run("Example restart-pa/SRL2: Serial-iter",log,value)

log = "./srl2Log/pa.err.1"
value = [['dragy',5.09547531705E-02,1e-06,4]]
Run("Example restart-pa/SRL2: Serial-error",log,value)



print("\n\pb:::")
#MPI
if ifmpi:

    log = "./mpiLog/pb.log.1"
    value = [['gmres: ',0,77,7]]
    Run("Example restart-pb/MPI: Serial-iter",log,value)

    log = "./mpiLog/pb.err.1"
    value = [['dragy',5.37986119139E-03,1e-06,4]]
    Run("Example restart-pb/MPI: Serial-error",log,value)

    log = "./mpiLog/pb.log.4"
    value = [['gmres: ',0,77,7]]
    Run("Example restart-pb/MPI: Parallel-iter",log,value)

    log = "./mpiLog/pb.err.4"
    value = [['dragy',5.37986119139E-03,1e-06,4]]
    Run("Example restart-pb/MPI: Parallel-error",log,value)

#SRL
log = "./srlLog/pb.log.1"
value = [['gmres: ',0,77,7]]
Run("Example restart-pb/SRL: Serial-iter",log,value)

log = "./srlLog/pb.err.1"
value = [['dragy',5.37986119139E-03,1e-06,4]]
Run("Example restart-pb/SRL: Serial-error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/pb.log.1"
    value = [['gmres: ',0,28,6]]
    Run("Example restart-pb/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/pb.err.1"
    value = [['dragy',5.09547531705E-02,1e-06,4]]
    Run("Example restart-pb/MPI2: Serial-error",log,value)

    log = "./mpi2Log/pb.log.4"
    value = [['gmres: ',0,28,6]]
    Run("Example restart-pb/MPI2: Parallel-iter",log,value)

    log = "./mpi2Log/pb.err.4"
    value = [['dragy',5.09547531705E-02,1e-06,4]]
    Run("Example restart-pb/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/pb.log.1"
value = [['gmres: ',0,28,6]]
Run("Example restart-pb/SRL2: Serial-iter",log,value)

log = "./srl2Log/pb.err.1"
value = [['dragy',5.09547531705E-02,1e-06,4]]
Run("Example restart-pb/SRL2: Serial-error",log,value)




print("\n\neddy Example")  
#MPI
if ifmpi:

    log = "./mpiLog/eddy_uv.log.1"
    value = [['gmres: ',0,34,7]]
    Run("Example eddy/MPI: Serial-iter",log,value)

    log = "./mpiLog/eddy_uv.err.1"
    value = [['X err',6.007702E-07,1e-06,6],['Y err',6.489061E-07,1e-06,6]]
    Run("Example eddy/MPI: Serial-error",log,value)

    log = "./mpiLog/eddy_uv.log.4"
    value = [['gmres: ',0,34,7]]
    Run("Example eddy/MPI: Parallel-iter",log,value)

    log = "./mpiLog/eddy_uv.err.4"
    value = [['X err',6.007702E-07,1e-06,6],['Y err',6.489061E-07,1e-06,6]]
    Run("Example eddy/MPI: Parallel-error",log,value)

#SRL
log = "./srlLog/eddy_uv.log.1"
value = [['total solver time',0.1,80,2],
         ['gmres: ',0,34,7]]
Run("Example eddy/SRL: Serial-time/iter",log,value)

log = "./srlLog/eddy_uv.err.1"
value = [['X err',6.007702E-07,1e-06,6],['Y err',6.489061E-07,1e-06,6]]
Run("Example eddy/SRL: Serial-error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/eddy_uv.log.1"
    value = [['gmres: ',0,22,6]]
    Run("Example eddy/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/eddy_uv.err.1"
    value = [['X err',6.759103E-05,1e-06,6],['Y err',7.842019E-05,1e-06,6]]
    Run("Example eddy/MPI2: Serial-error",log,value)

    log = "./mpi2Log/eddy_uv.log.4"
    value = [['gmres: ',0,22,6]]
    Run("Example eddy/MPI2: Parallel-iter",log,value)

    log = "./mpi2Log/eddy_uv.err.4"
    value = [['X err',6.759103E-05,1e-06,6],['Y err',7.842019E-05,1e-06,6]]
    Run("Example eddy/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/eddy_uv.log.1"
value = [['total solver time',0.1,80,2],
         ['gmres: ',0,22,6]]
Run("Example eddy/SRL2: Serial-time/iter",log,value)

log = "./srl2Log/eddy_uv.err.1"
value = [['X err',6.759103E-05,1e-06,6],['Y err',7.842019E-05,1e-06,6]]
Run("Example eddy/SRL2: Serial-error",log,value)




print("\n\nAMG_eddy Example")  
#MPI
if ifmpi:

    log = "./mpiLog/amg_eddy.log.1"
    value = [['gmres: ',0,37,7]]
    Run("Example AMG_eddy/MPI: Serial-iter",log,value)

    log = "./mpiLog/amg_eddy.err.1"
    value = [['X err',6.007702E-07,1e-06,6],['Y err',6.489061E-07,1e-06,6]]
    Run("Example AMG_eddy/MPI: Serial-error",log,value)

    log = "./mpiLog/amg_eddy.log.4"
    value = [['gmres: ',0,37,7]]
    Run("Example AMG_eddy/MPI: Parallel-iter",log,value)

    log = "./mpiLog/amg_eddy.err.4"
    value = [['X err',6.007702E-07,1e-06,6],['Y err',6.489061E-07,1e-06,6]]
    Run("Example AMG_eddy/MPI: Parallel-error",log,value)

#SRL
log = "./srlLog/amg_eddy.log.1"
value = [['total solver time',0.1,80,2],
         ['gmres: ',0,37,7]]
Run("Example AMG_eddy/SRL: Serial-time/iter",log,value)

log = "./srlLog/amg_eddy.err.1"
value = [['X err',6.007702E-07,1e-06,6],['Y err',6.489061E-07,1e-06,6]]
Run("Example AMG_eddy/SRL: Serial-error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/amg_eddy.log.1"
    value = [['gmres: ',0,37,6]]
    Run("Example AMG_eddy/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/amg_eddy.err.1"
    value = [['X err',6.759103E-05,1e-06,6],['Y err',7.842019E-05,1e-06,6]]
    Run("Example AMG_eddy/MPI2: Serial-error",log,value)

    log = "./mpi2Log/amg_eddy.log.4"
    value = [['gmres: ',0,37,6]]
    Run("Example AMG_eddy/MPI2: Parallel-iter",log,value)

    log = "./mpi2Log/amg_eddy.err.4"
    value = [['X err',6.759103E-05,1e-06,6],['Y err',7.842019E-05,1e-06,6]]
    Run("Example AMG_eddy/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/amg_eddy.log.1"
value = [['total solver time',0.1,120,2],
         ['gmres: ',0,37,6]]
Run("Example AMG_eddy/SRL2: Serial-time/iter",log,value)

log = "./srl2Log/amg_eddy.err.1"
value = [['X err',6.759103E-05,1e-06,6],['Y err',7.842019E-05,1e-06,6]]
Run("Example AMG_eddy/SRL2: Serial-error",log,value)



 
print("\n\nhpts_ed Example")  
#MPI
if ifmpi:

    log = "./mpiLog/hpts_ed.log.1"
    value = "ABORT: Too many pts to read in hpts"
    FindPhrase("Example hpts_ed/MPI: Serial",log,value)

    log = "./mpiLog/hpts_ed.log.4"
    value = [['gmres: ',0,34,7]]
    Run("Example hpts_ed/MPI: Parallel-iter",log,value)

    log = "./mpiLog/hpts_ed.err.4"
    value = [['X err',6.007702E-07,1e-06,6],['Y err',6.489061E-07,1e-06,6]]
    Run("Example hpts_ed/MPI: Parallel-error",log,value)

    log = "./mpiLog/hpts_ed.err.0"
    value = [['  1.0000000E-01  ',-4.396170E-01,1e-05,1]]
    Run("Example hpts_ed/MPI: Parallel-HPTS",log,value)

#SRL
log = "./srlLog/hpts_ed.log.1"
value = [['total solver time',0.1,80,2],
         ['gmres: ',0,34,7]]
Run("Example hpts_ed/SRL: Serial-time/iter",log,value)

log = "./srlLog/hpts_ed.err.1"
value = [['X err',6.007702E-07,1e-06,6],['Y err',6.489061E-07,1e-06,6]]
Run("Example hpts_ed/SRL: Serial-error",log,value)

log = "./srlLog/hpts_ed.err.0"
value = [['  1.0000000E-01  ',-4.396170E-01,1e-05,1]]
Run("Example hpts_ed/SRL: SRL-HPTS",log,value)


#MPI2
if ifmpi:

    log = "./mpi2Log/hpts_ed.log.1"
    value = "ABORT: Too many pts to read in hpts"
    FindPhrase("Example hpts_ed/MPI2: Serial",log,value)

    log = "./mpi2Log/hpts_ed.log.4"
    value = [['gmres: ',0,22,6]]
    Run("Example hpts_ed/MPI2: Parallel-iter",log,value)

    log = "./mpi2Log/hpts_ed.err.4"
    value = [['X err',6.759103E-05,1e-06,6],['Y err',7.842019E-05,1e-06,6]]
    Run("Example hpts_ed/MPI2: Parallel-error",log,value)

    log = "./mpi2Log/hpts_ed.err.0"
    value = [['  1.0000000E-01  ',-4.3953660E-01,1e-06,1]]
    Run("Example hpts_ed/MPI2: Parallel-HPTS",log,value)

#SRL2
log = "./srl2Log/hpts_ed.log.1"
value = [['total solver time',0.1,80,2],
         ['gmres: ',0,22,6]]
Run("Example hpts_ed/SRL2: Serial-time/iter",log,value)

log = "./srl2Log/hpts_ed.err.1"
value = [['X err',6.759103E-05,1e-06,6],['Y err',7.842019E-05,1e-06,6]]
Run("Example hpts_ed/SRL2: Serial-error",log,value)

log = "./srl2Log/hpts_ed.err.0"
value = [['  1.0000000E-01  ',-4.3953660E-01,1e-06,1]]
Run("Example hpts_ed/SRL2: SRL-HPTS",log,value)




print("\n\neddy_neknek Example")  
#MPI
if ifmpi:

    log = "./mpiLog/eddy_neknek.err.2"
    value = [['X err   global' ,4.595605E-04,1e-06,7],
             ['Y err   global' ,6.887475E-04,1e-06,7]]
    Run("Example eddy_neknek/MPI: 2--error",log,value)

    log = "./mpiLog/eddy_neknek.err.4"
    value = [['X err   global' ,4.595605E-04,1e-06,7],
             ['Y err   global' ,6.887475E-04,1e-06,7]]
    Run("Example eddy_neknek/MPI: 4--error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/eddy_neknek.err.2"
    value = [['X err   global' ,3.925388E-03,1e-06,7],
             ['Y err   global',6.299443E-03,1e-06,7]]
    Run("Example eddy_neknek/MPI2: 2--error",log,value)

    log = "./mpi2Log/eddy_neknek.err.4"
    value = [['X err   global' ,3.925388E-03,1e-06,7],
             ['Y err   global',6.299443E-03,1e-06,7]]
    Run("Example eddy_neknek/MPI2: 4--error",log,value)




print("\n\neddy_psi_omega Example")  
#MPI
if ifmpi:

    log = "./mpiLog/psi_omega.err.1"
    value = [['X err' ,1.177007E-10,1e-06,6]]
    Run("Example Eddy psi_omega/MPI: 2--error",log,value)

    log = "./mpiLog/psi_omega.err.4"
    value = [['X err' ,1.177007E-10,1e-06,6]]
    Run("Example Eddy psi_omega/MPI: 4--error",log,value)

#SRL
log = "./srlLog/psi_omega.log.1"
value = [['total solver time',0.1,17,2]]
Run("Example Eddy psi_omega/SRL: Serial-time",log,value)

log = "./srlLog/psi_omega.err.1"
value = [['X err',1.177007E-10,1e-06,6]]
Run("Example Eddy psi_omega/SRL: Serial-error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/psi_omega.err.1"
    value = [['X err' ,1.177007E-10,1e-06,6]]
    Run("Example Eddy psi_omega/MPI2: 1--error",log,value)

    log = "./mpi2Log/psi_omega.err.4"
    value = [['X err' ,1.177007E-10,1e-06,6]]
    Run("Example Eddy psi_omega/MPI2: 4--error",log,value)
 
#SRL2
log = "./srl2Log/psi_omega.log.1"
value = [['total solver time',0.1,17,2]]
Run("Example Eddy psi_omega/SRL2: Serial-time",log,value)

log = "./srl2Log/psi_omega.err.1"
value = [['X err',1.177007E-10,1e-06,6]]
Run("Example Eddy psi_omega/SRL2: Serial-error",log,value)




print("\n\nexpansion Example")  
#MPI 
if ifmpi:

    log = "./mpiLog/expansion.log.1"
    value = [['gmres: ',0,77,7]]
    Run("Example expansion/MPI: Serial-iter",log,value)

    log = "./mpiLog/expansion.err.1"
    value = [['ubar',2.0087E+00,1e-06,2]]
    Run("Example expansion/MPI: Serial-error",log,value)

    log = "./mpiLog/expansion.log.4"
    value = [['gmres: ',0,77,7]]
    Run("Example expansion/MPI: Serial-iter",log,value)

    log = "./mpiLog/expansion.err.4"
    value = [['ubar',2.0087E+00,1e-06,2]]
    Run("Example expansion/MPI: Serial-error",log,value)

#SRL
log = "./srlLog/expansion.log.1"
value = [['total solver time',0.1,250,2],
         ['gmres: ',0,77,7]]
Run("Example expansion/SRL: Serial-time/iter",log,value)

log = "./srlLog/expansion.err.1"
value = [['ubar',2.0087E+00,1e-06,2]]
Run("Example ext_cyl/SRL: Serial-error",log,value)

#MPI2 
if ifmpi:

    log = "./mpi2Log/expansion.log.1"
    value = [['gmres: ',0,70,6]]
    Run("Example expansion/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/expansion.err.1"
    value = [['ubar',2.0000E+00,1e-06,2]]
    Run("Example expansion/MPI2: Serial-error",log,value)

    log = "./mpi2Log/expansion.log.4"
    value = [['gmres: ',0,70,6]]
    Run("Example expansion/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/expansion.err.4"
    value = [['ubar',2.0000E+00,1e-06,2]]
    Run("Example expansion/MPI2: Serial-error",log,value)

#SRL
log = "./srl2Log/expansion.log.1"
value = [['total solver time',0.1,150,2],
         ['gmres: ',0,70,6]]
Run("Example expansion/SRL2: Serial-time/iter",log,value)

log = "./srl2Log/expansion.err.1"
value = [['ubar',2.0000E+00,1e-06,2]]
Run("Example ext_cyl/SRL2: Serial-error",log,value)




print("\n\next_cyl Example")  
#MPI 
if ifmpi:

    log = "./mpiLog/ext_cyl.log.1"
    value = [['gmres: ',0,85,7]]
    Run("Example ext_cyl/MPI: Serial-iter",log,value)

    log = "./mpiLog/ext_cyl.err.1"
    value = [['dragx',1.2138790E+00,1e-06,4],['dragy',1.3040301E-07,1e-06,4]]
    Run("Example ext_cyl/MPI: Serial-error",log,value)

    log = "./mpiLog/ext_cyl.log.4"
    value = [['gmres: ',0,85,7]]
    Run("Example ext_cyl/MPI: Parallel-iter",log,value)

    log = "./mpiLog/ext_cyl.err.4"
    value = [['dragx',1.2138790E+00,1e-06,4],['dragy',1.3040301E-07,1e-06,4]]
    Run("Example ext_cyl/MPI: Parallel-error",log,value)

#SRL
log = "./srlLog/ext_cyl.log.1"
value = [['total solver time',0.1,400,2],
         ['gmres: ',0,85,7]]
Run("Example ext_cyl/SRL: Serial-time/iter",log,value)

log = "./srlLog/ext_cyl.err.1"
value = [['dragx',1.2138790E+00,1e-06,4],['dragy',1.3040301E-07,1e-06,4]]
Run("Example ext_cyl/SRL: Serial-error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/ext_cyl.log.1"
    value = [['gmres: ',0,26,6]]
    Run("Example ext_cyl/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/ext_cyl.err.1"
    value = [['dragx',1.2138878E+00,1e-05,4],['dragy',3.2334222E-07,1e-06,4]]
    Run("Example ext_cyl/MPI2: Serial-error",log,value)

    log = "./mpi2Log/ext_cyl.log.4"
    value = [['gmres: ',0,26,6]]
    Run("Example ext_cyl/MPI2: Parallel-iter",log,value)

    log = "./mpi2Log/ext_cyl.err.4"
    value = [['dragx',1.2138878E+00,1e-05,4],['dragy',3.2334222E-07,1e-06,4]]
    Run("Example ext_cyl/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/ext_cyl.log.1"
value = [['total solver time',0.1,380,2],
         ['gmres: ',0,26,6]]
Run("Example ext_cyl/SRL2: Serial-time/iter",log,value)

log = "./srl2Log/ext_cyl.err.1"
value = [['dragx',1.2138878E+00,1e-05,4],['dragy',3.2334222E-07,1e-06,4]]
Run("Example ext_cyl/SRL2: Serial-error",log,value)




print("\n\nfs_2-st1 Example")  
#MPI 
if ifmpi:

    log = "./mpiLog/st1.log.1"
    value = "ABORT: "
    FindPhrase("Example st1/MPI: Serial",log,value)

    log = "./mpiLog/st1.log.4"
    value = "ABORT: "
    FindPhrase("Example st1/MPI: Parallel",log,value)

#SRL
log = "./srlLog/st1.log.1"
value = "ABORT: "
FindPhrase("Example st1/SRL: Serial",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/st1.log.1"
    value = [['gmres: ',0,38,6]]
    Run("Example st1/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/st1.err.1"
    value = [['amp',6.382414E-01,1e-06,2]]
    Run("Example st1/MPI2: Serial-error",log,value)

    log = "./mpi2Log/st1.log.4"
    value = [['gmres: ',0,38,6]]
    Run("Example st1/MPI2: Parallel-iter",log,value)

    log = "./mpi2Log/st1.err.4"
    value = [['amp',6.382414E-01,1e-06,2]]
    Run("Example st1/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/st1.log.1"
value = [['total solver time',0.1,18.3,2],
         ['gmres: ',0,38,6]]
Run("Example st1/SRL2: Serial-time/iter",log,value)

log = "./srl2Log/st1.err.1"
value = [['amp',6.382414E-01,1e-06,2]]
Run("Example st1/SRL2: Serial-error",log,value)




print("\n\nfs_2-st2 Example")  
#MPI 
if ifmpi:

    log = "./mpiLog/st2.log.1"
    value = "ABORT: "
    FindPhrase("Example st2/MPI: Serial",log,value)

    log = "./mpiLog/st2.log.4"
    value = "ABORT: "
    FindPhrase("Example st2/MPI: Parallel",log,value)

#SRL
log = "./srlLog/st2.log.1"
value = "ABORT: "
FindPhrase("Example st2/SRL: Serial",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/st2.log.1"
    value = [['gmres: ',0,38,6]]
    Run("Example st2/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/st2.err.1"
    value = [['amp',6.376171E-01,1e-06,2]]
    Run("Example st2/MPI2: Serial-error",log,value)

    log = "./mpi2Log/st2.log.4"
    value = [['gmres: ',0,38,6]]
    Run("Example st2/MPI2: Parallel-iter",log,value)

    log = "./mpi2Log/st2.err.4"
    value = [['amp',6.376171E-01,1e-06,2]]
    Run("Example st2/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/st2.log.1"
value = [['total solver time',0.1,23,2],
         ['gmres: ',0,38,6]]
Run("Example st2/SRL2: Serial-time/iter",log,value)

log = "./srl2Log/st2.err.1"
value = [['amp',6.376171E-01,1e-06,2]]
Run("Example st2/SRL2: Serial-error",log,value)




print("\n\nfs_2-std_wv Example")  
#MPI 
if ifmpi:

    log = "./mpiLog/std_wv.log.1"
    value = "ABORT: "
    FindPhrase("Example std_wv/MPI: Serial",log,value)

    log = "./mpiLog/std_wv.log.4"
    value = "ABORT: "
    FindPhrase("Example std_wv/MPI: Parallel",log,value)

#SRL
log = "./srlLog/std_wv.log.1"
value = "ABORT: "
FindPhrase("Example std_wv/SRL: Serial",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/std_wv.log.1"
    value = [['gmres: ',0,20,6]]
    Run("Example std_wv/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/std_wv.err.1"
    value = [['amp',1.403287E-01,1e-06,2]]
    Run("Example std_wv/MPI2: Serial-error",log,value)

    log = "./mpi2Log/std_wv.log.4"
    value = [['gmres: ',0,20,6]]
    Run("Example std_wv/MPI2: Parallel-iter",log,value)

    log = "./mpi2Log/std_wv.err.4"
    value = [['amp',1.403287E-01,1e-06,2]]
    Run("Example std_wv/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/std_wv.log.1"
value = [['total solver time',0.1,21,2],
         ['gmres: ',0,20,6]]
Run("Example std_wv/SRL2: Serial-time/iter",log,value)

log = "./srl2Log/std_wv.err.1"
value = [['amp',1.403287E-01,1e-06,2]]
Run("Example std_wv/SRL2: Serial-error",log,value)




print("\n\nfs_hydro Example")  
#MPI
if ifmpi:

    log = "./mpiLog/fs_hydro.log.1"
    value = "ABORT: "
    FindPhrase("Example fs_hydro/MPI: Serial",log,value)

    log = "./mpiLog/fs_hydro.log.4"
    value = "ABORT: "
    FindPhrase("Example fs_hydro/MPI: Parallel",log,value)

#SRL
log = "./srlLog/fs_hydro.log.1"
value = "ABORT: "
FindPhrase("Example fs_hydro/SRL: Serial",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/fs_hydro.log.1"
    value = [['gmres: ',0,108,6]]
    Run("Example fs_hydro/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/fs_hydro.err.1"
    value = [['AMP',-6.4616452E-05,2e-03,2]]
    Run("Example fs_hydro/MPI2: Serial-error",log,value)

    log = "./mpi2Log/fs_hydro.log.4"
    value = [['gmres: ',0,108,6]]
    Run("Example fs_hydro/MPI2: Parallel-iter",log,value)

    log = "./mpi2Log/fs_hydro.err.4"
    value = [['AMP',-6.4616452E-05,2e-03,2]]
    Run("Example fs_hydro/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/fs_hydro.log.1"
value = [['total solver time',0.1,200,2],
         ['gmres: ',0,108,6]]
Run("Example fs_hydro/SRL2: Serial-time/iter",log,value)

log = "./srl2Log/fs_hydro.err.1"
value = [['AMP',-6.4616452E-05,2e-03,2]]
Run("Example fs_hydro/SRL2: Serial-error",log,value)




print("\n\nhemi Example")  
#MPI
if ifmpi:

    log = "./mpiLog/hemi.log.1"
    value = [['gmres: ',0,39,7]]
    Run("Example hemi/MPI: Serial-iter",log,value)

    log = "./mpiLog/hemi.err.1"
    value = [['wmax',4.9173E-01,1e-06,2]]
    Run("Example hemi/MPI: Serial-error",log,value)

    log = "./mpiLog/hemi.log.4"
    value = [['gmres: ',0,39,7]]
    Run("Example hemi/MPI: Parallel-iter",log,value)

    log = "./mpiLog/hemi.err.4"
    value = [['wmax',4.9173E-01,1e-06,2]]
    Run("Example hemi/MPI: Parallel-error",log,value)

#SRL
log = "./srlLog/hemi.log.1"
value = [['total solver time',0.1,100,2],
         ['gmres: ',0,39,7]]
Run("Example hemi/SRL: Serial-time/iter",log,value)

log = "./srlLog/hemi.err.1"
value = [['wmax',4.9173E-01,1e-06,2]]
Run("Example hemi/SRL: Serial-error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/hemi.log.1"
    value = [['gmres: ',0,34,6]]
    Run("Example hemi/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/hemi.err.1"
    value = [['wmax',4.7915E-01,1e-06,2]]
    Run("Example hemi/MPI2: Serial-error",log,value)

    log = "./mpi2Log/hemi.log.4"
    value = [['gmres: ',0,34,6]]
    Run("Example hemi/MPI2: Parallel-iter",log,value)

    log = "./mpi2Log/hemi.err.4"
    value = [['wmax',4.7915E-01,1e-06,2]]
    Run("Example hemi/MPI2: Parallel-error",log,value)

#SRL
log = "./srl2Log/hemi.log.1"
value = [['total solver time',0.1,60,2],
         ['gmres: ',0,34,6]]
Run("Example hemi/SRL2: Serial-time/iter",log,value)

log = "./srl2Log/hemi.err.1"
value = [['wmax',4.7915E-01,1e-06,2]]
Run("Example hemi/SRL2: Serial-error",log,value)




print("\n\nkovasznay Example")  
#MPI
if ifmpi:

    log = "./mpiLog/kov.log.1"
    value = [['gmres: ',0,34,7]]
    Run("Example kov/MPI: Serial-iter",log,value)

    log = "./mpiLog/kov.err.1"
    value = [['err',5.14316E-13,1e-06,3]]
    Run("Example kov/MPI: Serial-error",log,value)

    log = "./mpiLog/kov.log.4"
    value = [['gmres: ',0,34,7]]
    Run("Example kov/MPI: Parallel-iter",log,value)

    log = "./mpiLog/kov.err.4"
    value = [['err',5.14316E-13,1e-06,3]]
    Run("Example kov/MPI: Parallel-error",log,value)

#SRL
log = "./srlLog/kov.log.1"
value = [['total solver time',0.1,12,2],
         ['gmres: ',0,34,7]]
Run("Example kov/SRL: Serial-time/iter",log,value)

log = "./srlLog/kov.err.1"
value = [['err',5.14316E-13,1e-06,3]]
Run("Example kov/SRL: Serial-error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/kov.log.1"
    value = [['gmres: ',0,14,6]]
    Run("Example kov/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/kov.err.1"
    value = [['err',5.90551E-13,1e-06,3]]
    Run("Example kov/MPI2: Serial-error",log,value)

    log = "./mpi2Log/kov.log.4"
    value = [['gmres: ',0,14,6]]
    Run("Example kov/MPI2: Parallel-iter",log,value)

    log = "./mpi2Log/kov.err.4"
    value = [['err',5.90551E-13,1e-06,3]]
    Run("Example kov/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/kov.log.1"
value = [['total solver time',0.1,17,2],
         ['gmres: ',0,14,6]]
Run("Example kov/SRL2: Serial-time/iter",log,value)

log = "./srl2Log/kov.err.1"
value = [['err',5.90551E-13,1e-06,3]]
Run("Example kov/SRL2: Serial-error",log,value)




print("\n\nkov_st_state Example")  
#MPI2
if ifmpi:

    log = "./mpi2Log/kov_st_stokes.err.1"
    value = [['err',8.55641E-10,1e-06,3]]
    Run("Example kov_st_state/MPI2: Serial-error",log,value)

    log = "./mpi2Log/kov_st_stokes.err.4"
    value = [['err',8.55641E-10,1e-06,3]]
    Run("Example kov_st_state/MPI2: Parallel-error",log,value)
 
#SRL2
log = "./srl2Log/kov_st_stokes.log.1"
value = [['total solver time',0.1,5,2]]
Run("Example kov_st_state/SRL2: Serial-time",log,value)

log = "./srl2Log/kov_st_stokes.err.1"
value = [['err',8.55641E-10,1e-06,3]]
Run("Example kov/SRL2: Serial-error",log,value)




print("\n\nlowMach_test Example")  
#MPI
if ifmpi:

    log = "./mpiLog/lowMach_test.log.1"
    value = [['gmres: ',0,100,7]]
    Run("Example lowMach_test/MPI: Serial-iter",log,value)

    log = "./mpiLog/lowMach_test.err.1"
    value = [['VX',2.4635E-09,1e-06,5],['T',4.5408E-12,1e-06,5],['QTL',2.6557E-06,1e-06,5]]
    Run("Example lowMach_test/MPI: Serial-error",log,value)

    log = "./mpiLog/lowMach_test.log.4"
    value = [['gmres: ',0,100,7]]
    Run("Example lowMach_test/MPI: Parallel-iter",log,value)

    log = "./mpiLog/lowMach_test.err.4"
    value = [['VX',2.4635E-09,1e-06,5],['T',4.5408E-12,1e-06,5],['QTL',2.6557E-06,1e-06,5]]
    Run("Example lowMach_test/MPI: Parallel-error",log,value)

#SRL
log = "./srlLog/lowMach_test.log.1"
value = [['total solver time',0.1,40,2],
         ['gmres: ',0,100,7]]
Run("Example lowMach_test/SRL: Serial-time/iter",log,value)

log = "./srlLog/lowMach_test.err.1"
value = [['VX',2.4635E-09,1e-06,5],['T',4.5408E-12,1e-06,5],['QTL',2.6557E-06,1e-06,5]]
Run("Example lowMach_test/SRL: Serial-error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/lowMach_test.log.1"
    value = "ABORT: For lowMach,"
    FindPhrase("Example lowMach_test/MPI2: Serial",log,value)

    log = "./mpi2Log/lowMach_test.log.4"
    value = "ABORT: For lowMach,"
    FindPhrase("Example lowMach_test/MPI2: Parallel",log,value)

#SRL2
log = "./srl2Log/lowMach_test.log.1"
value = "ABORT: For lowMach,"
FindPhrase("Example lowMach_test/SRL2: Serial",log,value)




print("\n\nmhd-gpf Example")  
#MPI
if ifmpi:

    log = "./mpiLog/gpf.log.1"
    value = "ABORT: MHD"
    FindPhrase("Example MHD-gpf/MPI: Serial",log,value)

    log = "./mpiLog/gpf.log.4"
    value = "ABORT: MHD"
    FindPhrase("Example MHD-gpf/MPI: Parallel",log,value)

#SRL
log = "./srlLog/gpf.log.1"
value = "ABORT: MHD"
FindPhrase("Example MHD-gpf/SRL: Serial",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/gpf.log.1"
    value = [['gmres: ',0,15,6]]
    Run("Example MHD-gpf/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/gpf.err.1"
    value = [['rtavg_gr_Em', 2.56712250E-01,.02,4]]
    Run("Example MHD-gpf/MPI2: Serial-error",log,value)

    log = "./mpi2Log/gpf.log.4"
    value = [['gmres: ',0,15,6]]
    Run("Example MHD-gpf/MPI2: Parallel-iter",log,value)

    log = "./mpi2Log/gpf.err.4"
    value = [['rtavg_gr_Em', 2.56712250E-01,.02,4]]
    Run("Example MHD-gpf/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/gpf.log.1"
value = [['total solver time',0.1,130,2],
         ['gmres: ',0,15,6]]
Run("Example MHD-gpf/SRL2: Serial-time/iter",log,value)

log = "./srl2Log/gpf.err.1"
value = [['rtavg_gr_Em', 2.56712250E-01,.02,4]]
Run("Example MHD-gpf/SRL2: Serial-error",log,value)




print("\n\nmhd-gpf_m Example")  
#MPI
if ifmpi:

    log = "./mpiLog/gpf_m.log.1"
    value = "ERROR: FDM"
    FindPhrase("Example MHD-gpf_m/MPI: Serial",log,value)

    log = "./mpiLog/gpf_m.log.4"
    value = "ERROR: FDM"
    FindPhrase("Example MHD-gpf_m/MPI: Parallel",log,value)

#SRL
log = "./srlLog/gpf_m.log.1"
value = "ERROR: FDM"
FindPhrase("Example MHD-gpf_m/SRL: Serial",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/gpf_m.err.1"
    value = [['rtavg_gr_Em', 2.56712250E-01,.02,4]]
    Run("Example MHD-gpf_m/MPI2: Serial-error",log,value)

    log = "./mpi2Log/gpf_m.err.4"
    value = [['rtavg_gr_Em', 2.56712250E-01,.02,4]]
    Run("Example MHD-gpf_m/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/gpf_m.log.1"
value = [['total solver time',0.1,130,2]]
Run("Example MHD-gpf_m/SRL2: Serial-time",log,value)

log = "./srl2Log/gpf_m.err.1"
value = [['rtavg_gr_Em', 2.56712250E-01,.02,4]]
Run("Example MHD-gpf_m/SRL2: Serial-error",log,value)




print("\n\nmhd-gpf_b Example")  
#MPI
if ifmpi:

    log = "./mpiLog/gpf_b.log.1"
    value = "ABORT: MHD"
    FindPhrase("Example MHD-gpf_b/MPI: Serial",log,value)

    log = "./mpiLog/gpf_b.log.4"
    value = "ABORT: MHD"
    FindPhrase("Example MHD-gpf_b/MPI: Parallel",log,value)

#SRL
log = "./srlLog/gpf_b.log.1"
value = "ABORT: MHD"
FindPhrase("Example MHD-gpf_b/SRL: Serial",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/gpf_b.err.1"
    value = [['rtavg_gr_Em', 2.56712250E-01,.02,4]]
    Run("Example MHD-gpf_b/MPI2: Serial-error",log,value)

    log = "./mpi2Log/gpf_b.err.4"
    value = [['rtavg_gr_Em', 2.56712250E-01,.02,4]]
    Run("Example MHD-gpf_b/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/gpf_b.err.1"
value = [['rtavg_gr_Em', 2.56712250E-01,.02,4]]
Run("Example MHD-gpf_b/SRL2: Serial-error",log,value)




#print("\n\nmoab Example")  
#MPI
#if ifmpi:
#
#    log = "./mpiLog/pipe.log.1"
#    value = [['gmres: ',0,17,7]]
#    Run("Example moab/MPI: Serial-iter",log,value)
#
#    log = "./mpiLog/pipe.log.4"
#    value = [['gmres: ',0,17,7]]
#    Run("Example moab/MPI: Parallel-iter",log,value)
#
#SRL
#log = "./srlLog/pipe.log.1"
#value = [['total solver time',0.1,180,2],
#         ['gmres: ',0,17,7]]
#Run("Example moab/SRL: Serial-time/iter",log,value)
#
##MPI2
#if ifmpi:
#
#    log = "./mpi2Log/pipe.log.1"
#    value = "ABORT: "
#    FindPhrase("Example moab/MPI2: Serial",log,value)
#
#    log = "./mpi2Log/pipe.log.4"
#    value = "ABORT: "
#    FindPhrase("Example moab/MPI2: Parallel",log,value)
# 
##SRL2
#log = "./srl2Log/pipe.log.1"
#value = "ABORT: "
#FindPhrase("Example moab/SRL2: Serial",log,value)
#
#
##print("\n\nmoab_conjht Example")  
#MPI
#if ifmpi:
#
#    log = "./mpiLog/moab_conjht.log.1"
#    value = [['gmres: ',0,13,6]]
#    Run("Example moab_conjht/MPI: Serial-iter",log,value)
#
#    log = "./mpiLog/moab_conjht.err.1"
#    value = [['tmax',3.13266E+00,1e-06,2]]
#    Run("Example moab_conjht/MPI: Serial-error",log,value)
#
#    log = "./mpiLog/moab_conjht.log.2"
#    value = [['gmres: ',0,13,6]]
#    Run("Example moab_conjht/MPI: Parallel-iter",log,value)
#
#    log = "./mpiLog/moab_conjht.err.2"
#    value = [['tmax',3.13266E+00,1e-06,2]]
#    Run("Example moab_conjht/MPI: Parallel-error",log,value)
#
#SRL
#log = "./srlLog/moab_conjht.log.1"
#value = [['total solver time',0.1,10,2],
#         ['gmres: ',0,35,6]]
#Run("Example moab_conjht/SRL: Serial-time/iter",log,value)
#
#log = "./srlLog/moab_conjht.err.1"
#value = [['tmax',1.59864E+00,1e-06,2]]
#Run("Example moab_conjht/SRL: Serial-error",log,value)

#MPI2
#if ifmpi:
#
#    log = "./mpi2Log/moab_conjht.log.1"
#    value = "ABORT: "
#    FindPhrase("Example moab_conjht/MPI2: Serial",log,value)
#
#    log = "./mpi2Log/moab_conjht.log.2"
#    value = "ABORT: "
#    FindPhrase("Example moab_conjht/MPI2: Parallel-iter",log,value)
#
#
#log = "./srl2Log/moab_conjht.log.1"
#value = "ABORT: "
#FindPhrase("Example moab_conjht/SRL2: Serial-time/iter",log,value)




print("\n\nos7000 Example")  
#MPI
if ifmpi:

    log = "./mpiLog/u3_t020_n13.log.1"
    value = [['gmres: ',0,43,7]]
    Run("Example os7000/MPI: Serial-iter",log,value)

    log = "./mpiLog/u3_t020_n13.err.1"
    value = [['egn',4.74494769e-05,1e-06,2]]
    Run("Example os7000/MPI: Serial-error",log,value)

    log = "./mpiLog/u3_t020_n13.log.4"
    value = [['gmres: ',0,43,7]]
    Run("Example os7000/MPI: Parallel-iter",log,value)

    log = "./mpiLog/u3_t020_n13.err.4"
    value = [['egn',4.74494769e-05,1e-06,2]]
    Run("Example os7000/MPI: Parallel-error",log,value)

#SRL
log = "./srlLog/u3_t020_n13.log.1"
value = [['total solver time',0.1,40,2],
         ['gmres: ',0,43,7]]
Run("Example os7000/SRL: Serial-time/iter",log,value)

log = "./srlLog/u3_t020_n13.err.1"
value = [['egn',4.74494769e-05,1e-06,2]]
Run("Example os7000/SRL: Serial-error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/u3_t020_n13.log.1"
    value = [['gmres: ',0,43,6]]
    Run("Example os7000/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/u3_t020_n13.err.1"
    value = [['egn',5.93471252E-05,1e-06,2]]
    Run("Example os7000/MPI2: Serial-error",log,value)
    
    log = "./mpi2Log/u3_t020_n13.log.4"
    value = [['gmres: ',0,43,6]]
    Run("Example os7000/MPI2: Parallel-iter",log,value)

    log = "./mpi2Log/u3_t020_n13.err.4"
    value = [['egn',5.93471252E-05,1e-06,2]]
    Run("Example os7000/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/u3_t020_n13.log.1"
value = [['total solver time',0.1,40,2],
         ['gmres: ',0,43,6]]
Run("Example os7000/SRL2: Serial-iter",log,value)

log = "./srl2Log/u3_t020_n13.err.1"
value = [['egn',5.93471252E-05,1e-06,2]]
Run("Example os7000/SRL2: Serial-error",log,value)




print("\n\nperis Example")  
#MPI
if ifmpi:

    log = "./mpiLog/peris.log.1"
    value = "ABORT: "
    FindPhrase("Example peris/MPI: Serial",log,value)

    log = "./mpiLog/peris.log.4"
    value = "ABORT: "
    FindPhrase("Example peris/MPI: Parallel",log,value)

#SRL
log = "./srlLog/peris.log.1"
value = "ABORT: "
FindPhrase("Example peris/SRL: Serial",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/peris.log.1"
    value = [['gmres: ',0,18,6]]
    Run("Example peris/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/peris.log.4"
    value = [['gmres: ',0,18,6]]
    Run("Example peris/MPI2: Parallel-iter",log,value)

#SRL2
log = "./srl2Log/peris.log.1"
value = [['total solver time',0.1,13,2],
         ['gmres: ',0,18,6]]
Run("Example peris/SRL2: Serial-time/iter",log,value)




print("\n\npipe-helix Example")  
#MPI 
if ifmpi:

    log = "./mpiLog/helix.log.1"
    value = [['gmres: ',0,61,7]]
    Run("Example helix/MPI: Serial-iter",log,value)

    log = "./mpiLog/helix.err.1"
    value = [['err2',1.9077617E+00,1e-06,2]]
    Run("Example helix/MPI: Serial-error",log,value)

    log = "./mpiLog/helix.log.4"
    value = [['gmres: ',0,61,7]]
    Run("Example helix/MPI: Parallel-iter",log,value)

    log = "./mpiLog/helix.err.4"
    value = [['err2',1.9077617E+00,1e-06,2]]
    Run("Example helix/MPI: Serial-error",log,value)

#SRL
log = "./srlLog/helix.log.1"
value = [['total solver time',0.1,22,2],
         ['gmres: ',0,61,7]]
Run("Example helix/SRL: Serial-time/iter",log,value)

log = "./srlLog/helix.err.1"
value = [['err2',1.9077617E+00,1e-06,2]]
Run("Example helix/SRL: Serial-error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/helix.log.1"
    value = [['gmres: ',0,123,6]]
    Run("Example helix/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/helix.err.1"
    value = [['err2',1.9072258E+00,1e-06,2]]
    Run("Example helix/MPI2: Serial-error",log,value)

    log = "./mpi2Log/helix.log.4"
    value = [['gmres: ',0,123,6]]
    Run("Example helix/MPI2: Parallel-iter",log,value)

    log = "./mpi2Log/helix.err.4"
    value = [['err2',1.9072258E+00,1e-06,2]]
    Run("Example helix/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/helix.log.1"
value = [['total solver time',0.1,22,2],
         ['gmres: ',0,123,6]]
Run("Example helix/SRL2: Serial-time/iter",log,value)

log = "./srl2Log/helix.err.1"
value = [['err2',1.9072258E+00,1e-06,2]]
Run("Example helix/SRL2: Serial-error",log,value)




print("\n\npipe-stenosis Example")  
#MPI 
if ifmpi:

    log = "./mpiLog/stenosis.log.1"
    value = [['gmres: ',0,196,7]]
    Run("Example stenosis/MPI: Serial-iter",log,value)

    log = "./mpiLog/stenosis.log.4"
    value = [['gmres: ',0,196,7]]
    Run("Example stenosis/MPI: Parallel-iter",log,value)

#SRL
log = "./srlLog/stenosis.log.1"
value = [['total solver time',0.1,80,2],
         ['gmres: ',0,196,7]]
Run("Example stenosis/SRL: Serial-time/iter",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/stenosis.log.1"
    value = [['gmres: ',0,51,6]]
    Run("Example stenosis/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/stenosis.log.4"
    value = [['gmres: ',0,51,6]]
    Run("Example stenosis/MPI2: Parallel-iter",log,value)

#SRL2
log = "./srl2Log/stenosis.log.1"
value = [['total solver time',0.1,40,2],
         ['gmres: ',0,51,6]]
Run("Example stenosis/SRL2: Serial-time/iter",log,value)




print("\n\nrayleigh-ray1 Example")  
#MPI
if ifmpi:

    log = "./mpiLog/ray1.log.1"
    value = [['gmres: ',0,32,7]]
    Run("Example ray1/MPI: Serial-iter",log,value)

    log = "./mpiLog/ray1.err.1"
    value = [['umax',2.792052E-03,1e-03,3]]
    Run("Example ray1/MPI: Serial-error",log,value)

    log = "./mpiLog/ray1.log.4"
    value = [['gmres: ',0,32,7]]
    Run("Example ray1/MPI: Parallel-iter",log,value)

    log = "./mpiLog/ray1.err.4"
    value = [['umax',2.792052E-03,1e-03,3]]
    Run("Example ray1/MPI: Parallel-error",log,value)

#SRL
log = "./srlLog/ray1.log.1"
value = [['total solver time',0.1,3,2],
         ['gmres: ',0,32,7]]
Run("Example ray1/SRL: Serial-time/iter",log,value)

log = "./srlLog/ray1.err.1"
value = [['umax',2.792052E-03,1e-03,3]]
Run("Example ray1/SRL: Serial-error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/ray1.log.1"
    value = [['gmres: ',0,11,6]]
    Run("Example ray1/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/ray1.err.1"
    value = [['umax',4.831113E-03,1e-05,3]]
    Run("Example ray1/MPI2: Serial-error",log,value)

    log = "./mpi2Log/ray1.log.4"
    value = [['gmres: ',0,11,6]]
    Run("Example ray1/MPI2: Parallel-iter",log,value)

    log = "./mpi2Log/ray1.err.4"
    value = [['umax',4.831113E-03,1e-05,3]]
    Run("Example ray1/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/ray1.log.1"
value = [['total solver time',0.1,3,2],
         ['gmres: ',0,11,6]]
Run("Example ray1/SRL2: Serial-time/iter",log,value)

log = "./srl2Log/ray1.err.1"
value = [['umax',4.831113E-03,1e-05,3]]
Run("Example ray1/SRL2: Serial-error",log,value)




print("\n\nrayleigh-ray2 Example")  
#MPI
if ifmpi:

    log = "./mpiLog/ray2.log.1"
    value = [['gmres: ',0,31,7]]
    Run("Example ray2/MPI: Serial-iter",log,value)

    log = "./mpiLog/ray2.err.1"
    value = [['umax',4.549071E-03,1e-05,3]]
    Run("Example ray2/MPI: Serial-error",log,value)

    log = "./mpiLog/ray2.log.4"
    value = [['gmres: ',0,31,7]]
    Run("Example ray2/MPI: Parallel-iter",log,value)

    log = "./mpiLog/ray2.err.4"
    value = [['umax',4.549071E-03,1e-05,3]]
    Run("Example ray2/MPI: Parallel-error",log,value)

#SRL
log = "./srlLog/ray2.log.1"
value = [['total solver time',0.1,3,2],
         ['gmres: ',0,31,7]]
Run("Example ray2/SRL: Serial-time/iter",log,value)

log = "./srlLog/ray2.err.1"
value = [['umax',4.549071E-03,1e-05,3]]
Run("Example ray2/SRL: Serial-error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/ray2.log.1"
    value = [['gmres: ',0,11,6]]
    Run("Example ray2/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/ray2.err.1"
    value = [['umax',6.728787E-03,1e-05,3]]
    Run("Example ray2/MPI2: Serial-error",log,value)

    log = "./mpi2Log/ray2.log.4"
    value = [['gmres: ',0,11,6]]
    Run("Example ray2/MPI2: Parellel-iter",log,value)

    log = "./mpi2Log/ray2.err.4"
    value = [['umax',6.728787E-03,1e-05,3]]
    Run("Example ray2/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/ray2.log.1"
value = [['total solver time',0.1,3,2],
         ['gmres: ',0,11,6]]
Run("Example ray2/SRL2: Serial-time/iter",log,value)

log = "./srl2Log/ray2.err.1"
value = [['umax',6.728787E-03,1e-05,3]]
Run("Example ray2/SRL2: Serial-error",log,value)




print("\n\nshear4-shear4 Example")  
#MPI 
if ifmpi:

    log = "./mpiLog/shear4.log.1"
    value = [['gmres: ',0,26,7]]
    Run("Example shear4/thick/MPI: Serial-iter",log,value)

    log = "./mpiLog/shear4.err.1"
    value = [['peak vorticity',3.031328E+01,1e-06,3]]
    Run("Example shear4/thick/MPI: Serial-error",log,value)

    log = "./mpiLog/shear4.log.4"
    value = [['gmres: ',0,26,7]]
    Run("Example shear4/thick/MPI: Parallel-iter",log,value)

    log = "./mpiLog/shear4.err.4"
    value = [['peak vorticity',3.031328E+01,1e-06,3]]
    Run("Example shear4/thick/MPI: Parallel-error",log,value)

#SRL
log = "./srlLog/shear4.log.1"
value = [['total solver time',0.1,10,2],
         ['gmres: ',0,26,7]]
Run("Example shear4/thick/SRL: Serial-time/iter",log,value)

log = "./srlLog/shear4.err.1"
value = [['peak vorticity',3.031328E+01,1e-06,3]]
Run("Example shear4/thick/SRL: Serial-error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/shear4.log.1"
    value = [['gmres: ',0,17,6]]
    Run("Example shear4/thick/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/shear4.err.1"
    value = [['peak vorticity',3.031328E+01,1e-06,3]]
    Run("Example shear4/thick/MPI2: Serial-error",log,value)

    log = "./mpi2Log/shear4.log.4"
    value = [['gmres: ',0,17,6]]
    Run("Example shear4/thick/MPI2: Parallel-iter",log,value)

    log = "./mpi2Log/shear4.err.4"
    value = [['peak vorticity',3.031328E+01,1e-06,3]]
    Run("Example shear4/thick/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/shear4.log.1"
value = [['total solver time',0.1,10,2],
         ['gmres: ',0,17,6]]
Run("Example shear4/thick/SRL2: Serial-time/iter",log,value)

log = "./srl2Log/shear4.err.1"
value = [['peak vorticity',3.031328E+01,1e-06,3]]
Run("Example shear4/thick/SRL2: Serial-error",log,value)




print("\n\nshear4-thin Example")  
#MPI
if ifmpi:

    log = "./mpiLog/thin.log.1"
    value = [['gmres: ',0,26,7]]
    Run("Example shear4/thin/MPI: Serial-iter",log,value)

    log = "./mpiLog/thin.err.1"
    value = [['peak vorticity',9.991753E+01,1e-06,3]]
    Run("Example shear4/thin/MPI: Serial-error",log,value)

    log = "./mpiLog/thin.log.4"
    value = [['gmres: ',0,26,7]]
    Run("Example shear4/thin/MPI: Parallel-iter",log,value)

    log = "./mpiLog/thin.err.4"
    value = [['peak vorticity',9.991753E+01,1e-06,3]]
    Run("Example shear4/thin/MPI: Parallel-error",log,value)

#SRL
log = "./srlLog/thin.log.1"
value = [['total solver time',0.1,10,2],
         ['gmres: ',0,26,7]]
Run("Example shear4/thin/SRL: Serial-time/iter",log,value)

log = "./srlLog/thin.err.1"
value = [['peak vorticity',9.991753E+01,1e-06,3]]
Run("Example shear4/thin/SRL: Serial-error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/thin.log.1"
    value = [['gmres: ',0,17,6]]
    Run("Example shear4/thin/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/thin.err.1"
    value = [['peak vorticity',9.991556E+01,1e-06,3]]
    Run("Example shear4/thin/MPI2: Serial-error",log,value)

    log = "./mpi2Log/thin.log.4"
    value = [['gmres: ',0,17,6]]
    Run("Example shear4/thin/MPI2: Parallel-iter",log,value)

    log = "./mpi2Log/thin.err.4"
    value = [['peak vorticity',9.991556E+01,1e-06,3]]
    Run("Example shear4/thin/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/thin.log.1"
value = [['total solver time',0.1,10,2],
         ['gmres: ',0,17,6]]
Run("Example shear4/thin/SRL2: Serial-time/iter",log,value)

log = "./srl2Log/thin.err.1"
value = [['peak vorticity',9.991556E+01,1e-06,3]]
Run("Example shear4/thin/SRL2: Serial-error",log,value)




print("\n\nSolid Example")  
#MPI
if ifmpi:

    log = "./mpiLog/solid.log.1"
    value = "ABORT: "
    FindPhrase("Example solid/MPI: Serial",log,value)

    log = "./mpiLog/solid.log.4"
    value = "ABORT: "
    FindPhrase("Example solid/MPI: Parallel",log,value)

#SRL
log = "./srlLog/solid.log.1"
value = "ABORT: "
FindPhrase("Example solid/SRL: Serial",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/solid.err.1"
    value = [['error',7.821228E-05,1e-06,2]]
    Run("Example solid/MPI2: Serial-error",log,value)

    log = "./mpi2Log/solid.err.4"
    value = [['error',7.821228E-05,1e-06,2]]
    Run("Example solid/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/solid.err.1"
value = [['error',7.821228E-05,1e-06,2]]
Run("Example solid/SRL2: Serial-error",log,value)




print("\n\nStrat Example")  
#MPI
if ifmpi:

    log = "./mpiLog/re10f1000p0001.log.1"
    value = [['gmres: ',0,60,7]]
    Run("Example strat-01/MPI: Serial-iter",log,value)

    log = "./mpiLog/re10f1000p0001.log.4"
    value = [['gmres: ',0,60,7]]
    Run("Example strat-01/MPI: PARALLEL-iter",log,value)

    log = "./mpiLog/re10f1000p1000.log.1"
    value = [['gmres: ',0,60,7]]
    Run("Example strat-1000/MPI: Serial-iter",log,value)

    log = "./mpiLog/re10f1000p1000.log.4"
    value = [['gmres: ',0,60,7]]
    Run("Example strat-1000/MPI: PARALLEL-iter",log,value)

#SRL
log = "./srlLog/re10f1000p0001.log.1"
value = [['total solver time',0.1,140,2],
         ['gmres: ',0,60,7]]
Run("Example strat-01/SRL: Serial-time/iter",log,value)

log = "./srlLog/re10f1000p1000.log.1"
value = [['total solver time',0.1,140,2],
         ['gmres: ',0,60,7]]
Run("Example strat-1000/SRL: Serial-time/iter",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/re10f1000p0001.log.1"
    value = [['U-PRES ',0,27,6]]
    Run("Example strat-01/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/re10f1000p0001.log.4"
    value = [['U-PRES ',0,27,6]]
    Run("Example strat-01/MPI2: PARALLEL-iter",log,value)

    log = "./mpi2Log/re10f1000p1000.log.1"
    value = [['U-PRES ',0,27,6]]
    Run("Example strat-1000/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/re10f1000p1000.log.4"
    value = [['U-PRES ',0,27,6]]
    Run("Example strat-1000/MPI2: PARALLEL-iter",log,value)

#SRL2
log = "./srl2Log/re10f1000p0001.log.1"
value = [['total solver time',0.1,80,2],
         ['U-PRES ',0,27,6]]
Run("Example strat-01/SRL: Serial-time/iter",log,value)

log = "./srl2Log/re10f1000p1000.log.1"
value = [['total solver time',0.1,80,2],
         ['U-PRES ',0,27,6]]
Run("Example strat-1000/SRL: Serial-time/iter",log,value)




print("\n\nTaylor Example")  
#MPI
if ifmpi:

    log = "./mpiLog/taylor.log.1"
    value = [['gmres: ',0,23,7]]
    Run("Example taylor/MPI: Serial-iter",log,value)

    log = "./mpiLog/taylor.err.1"
    value = [['tq',4.13037E-06,1e-06,5],
             ['err',2.973648E-09,1e-06,2]]
    Run("Example taylor/MPI: Serial-error",log,value)

    log = "./mpiLog/taylor.log.4"
    value = [['gmres: ',0,23,7]]
    Run("Example taylor/MPI: Parallel-iter",log,value)

    log = "./mpiLog/taylor.err.4"
    value = [['tq',4.13037E-06,1e-06,5],
             ['err',2.973648E-09,1e-06,2]]
    Run("Example taylor/MPI: Parallel-error",log,value)

#SRL
log = "./srlLog/taylor.log.1"
value = [['total solver time',0.1,40,2],
         ['gmres: ',0,23,7]]
Run("Example taylor/SRL: Serial-time/iter",log,value)

log = "./srlLog/taylor.err.1"
value = [['tq',4.13037E-06,1e-06,5],
         ['err',2.973648E-09,1e-06,2]]
Run("Example taylor/SRL: Serial-error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/taylor.log.1"
    value = [['gmres: ',0,14,6]]
    Run("Example taylor/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/taylor.err.1"
    value = [['tq',4.10783E-06,1e-06,5],
             ['err',2.826284E-10,1e-06,2]]
    Run("Example taylor/MPI2: Serial-error",log,value)

    log = "./mpi2Log/taylor.log.4"
    value = [['gmres: ',0,14,6]]
    Run("Example taylor/MPI2: Parallel-iter",log,value)

    log = "./mpi2Log/taylor.err.4"
    value = [['tq',4.10783E-06,1e-06,5],
             ['err',2.826284E-10,1e-06,2]]
    Run("Example taylor/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/taylor.log.1"
value = [['total solver time',0.1,40,2],
         ['gmres: ',0,14,6]]
Run("Example taylor/SRL2: Serial-time/iter",log,value)

log = "./srl2Log/taylor.err.1"
value = [['tq',4.10783E-06,1e-06,5],
         ['err',2.826284E-10,1e-06,2]]
Run("Example taylor/SRL2: Serial-error",log,value)




print("\n\nturbChannel Example")  
#MPI
if ifmpi:

    log = "./mpiLog/turbChannel.log.1"
    value = [['gmres: ',0,95,7]]
    Run("Example turbChannel/MPI: Serial-iter",log,value)

    log = "./mpiLog/turbChannel.log.4"
    value = [['gmres: ',0,95,7]]
    Run("Example turbChannel/MPI: Parallel-iter",log,value)

#SRL
log = "./srlLog/turbChannel.log.1"
value = [['total solver time',0.1,200,2],
         ['gmres: ',0,95,7]]
Run("Example turbChannel/SRL: Serial-time/iter",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/turbChannel.log.1"
    value = [['gmres: ',0,26,6]]
    Run("Example turbChannel/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/turbChannel.log.4"
    value = [['gmres: ',0,26,6]]
    Run("Example turbChannel/MPI2: Parallel-iter",log,value)

#SRL2
log = "./srl2Log/turbChannel.log.1"
value = [['total solver time',0.1,140,2],
         ['gmres: ',0,26,6]]
Run("Example turbChannel/SRL2: Serial-time/iter",log,value)




print("\n\nvar_vis Example")  
#MPI 
if ifmpi:

    log = "./mpiLog/var_vis.log.1"
    value = "ABORT: "
    FindPhrase("Example var_vis/MPI: Serial",log,value)

    log = "./mpiLog/var_vis.log.4"
    value = "ABORT: "
    FindPhrase("Example var_vis/MPI: Parallel",log,value)

#SRL
log = "./srlLog/var_vis.log.1"
value = "ABORT: "
FindPhrase("Example var_vis/SRL: Serial",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/var_vis.log.1"
    value = [['gmres: ',0,19,6]]
    Run("Example var_vis/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/var_vis.log.4"
    value = [['gmres: ',0,19,6]]
    Run("Example var_vis/MPI2: Parallel-iter",log,value)

#SRL2
log = "./srl2Log/var_vis.log.1"
value = [['total solver time',0.1,30,2],
         ['gmres: ',0,19,6]]
Run("Example var_vis/SRL2: Serial-time/iter",log,value)




print("\n\nvortex Example")  
#MPI
if ifmpi:

    log = "./mpiLog/r1854a.log.1"
    value = [['gmres: ',0,65,7]]
    Run("Example vortex/MPI: Serial-int",log,value)

    log = "./mpiLog/r1854a.err.1"
    value = [['VMIN',-1.910312E-03,1e-05,2]]
    Run("Example vortex/MPI: Serial-error",log,value)

    log = "./mpiLog/r1854a.log.4"
    value = [['gmres: ',0,65,7]]
    Run("Example vortex/MPI: Parallel-int",log,value)

    log = "./mpiLog/r1854a.err.4"
    value = [['VMIN',-1.910312E-03,1e-05,2]]
    Run("Example vortex/MPI: Parallel-error",log,value)

#SRL
log = "./srlLog/r1854a.log.1"
value = [['total solver time',0.1,60,2],
         ['gmres: ',0,65,7]]
Run("Example vortex/SRL: Serial-time/iter",log,value)

log = "./srlLog/r1854a.err.1"
value = [['VMIN',-1.910312E-03,1e-05,2]]
Run("Example vortex/SRL: Serial-error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/r1854a.log.1"
    value = [['gmres: ',0,18,6]]
    Run("Example vortex/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/r1854a.err.1"
    value = [['VMIN',-1.839120E-03,1e-06,2]]
    Run("Example vortex/MPI2: Serial-error",log,value)

    log = "./mpi2Log/r1854a.log.4"
    value = [['gmres: ',0,18,6]]
    Run("Example vortex/MPI2: Parallel-iter",log,value)

    log = "./mpi2Log/r1854a.err.4"
    value = [['VMIN',-1.839120E-03,1e-06,2]]
    Run("Example vortex/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/r1854a.log.1"
value = [['total solver time',0.1,50,2],
         ['gmres: ',0,18,6]]
Run("Example vortex/SRL2: Serial-time/iter",log,value)

log = "./srl2Log/r1854a.err.1"
value = [['VMIN',-1.839120E-03,1e-06,2]]
Run("Example vortex/SRL2: Serial-error",log,value)




print("\n\nvortex2 Example")  
#MPI
if ifmpi:
#first nine time steps fail in pressure

    log = "./mpiLog/v2d.log.1"
    value = [['PRES:  ',0,100,4]]
    Run("Example vortex2/MPI: Serial-iter",log,value)

    log = "./mpiLog/v2d.err.1"
    value = [['umin',-1.453402E-03,1e-03,2],
             ['torqx',-1.7399905E-07,1e-06,2]]
    Run("Example vortex2/MPI: Serial-error",log,value)

    log = "./mpiLog/v2d.log.4"
    value = [['PRES:  ',0,100,4]]
    Run("Example vortex2/MPI: Parallel-iter",log,value)

    log = "./mpiLog/v2d.err.4"
    value = [['umin',-1.453402E-03,1e-03,2],
             ['torqx',-1.7399905E-07,1e-06,2]]
    Run("Example vortex2/MPI: Parallel-error",log,value)

#SRL
log = "./srlLog/v2d.log.1"
value = [['total solver time',0.1,80,2],
         ['PRES: ',0,100,4]]
Run("Example vortex2/SRL: Serial-time/iter",log,value)

log = "./srlLog/v2d.err.1"
value = [['umin',-1.453402E-03,1e-03,2],
         ['torqx',-1.7399905E-07,1e-06,2]]
Run("Example vortex2/SRL: Serial-error",log,value)

#MPI2
if ifmpi:

    log = "./mpi2Log/v2d.log.1"
    value = [['U-Press ',0,100,5]]
    Run("Example vortex2/MPI2: Serial-iter",log,value)

    log = "./mpi2Log/v2d.err.1"
    value = [['umin',-2.448980E-03,1e-03,2],
             ['torqx',-1.6276138E-07,1e-06,2]]
    Run("Example vortex2/MPI2: Serial-error",log,value)

    log = "./mpi2Log/v2d.log.4"
    value = [['U-Press ',0,100,5]]
    Run("Example vortex2/MPI2: Parallel-iter",log,value)

    log = "./mpi2Log/v2d.err.4"
    value = [['umin',-2.448980E-03,1e-03,2],
             ['torqx',-1.6276138E-07,1e-06,2]]
    Run("Example vortex2/MPI2: Parallel-error",log,value)

#SRL2
log = "./srl2Log/v2d.log.1"
value = [['total solver time',0.1,80,2],
         ['U-Press ',0,100,5]]
Run("Example vortex2/SRL2: Serial-time/iter",log,value)

log = "./srl2Log/v2d.err.1"
value = [['umin',-2.448980E-03,1e-03,2],
         ['torqx',-1.6276138E-07,1e-06,2]]
Run("Example vortex2/SRL2: Serial-error",log,value)

###############################################################################
###############################################################################
print("\n\nTest Summary :     %i/%i tests were successful"%(num_success,num_test))
print("End of top-down testing")
######################################################################
