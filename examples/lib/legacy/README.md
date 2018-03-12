NekTests
========
Legacy Buildbot Tests
---------------------

### To Run Buildbot Manually On Any Platform(with at least 4 processors)
1.  Make sure you have a current set of the examples from the repo
    in your nek5_svn/examples directoy. (And your nek5_svn/tests directory
    is up to date as well)

2. `$ cd nek5_svn/tests`

3. `$ cp BB_RunTest your_RunTest`
   Edit your_RunTest with the compilers, MOAB directories, Matlab executable path
   that you want to test.  The current BB_RunTest is set to test 3 compilers, 
   one MPI implementation, along with the MOAB and AMG(matlab) tests.

   To run the Analysis step within the you_RunTest (as opposed to running it later, 
   separately, in step 5), include the commented out line:
  
      ./Analysis.py ${PERFORMED_TESTS} &> Analysis.log
  
   Remember to add a 
      mv Analysis.log ./$COMPILER
   to your script to move this log to the directoy of your test results.

   Note -- to disable any test, leave that parameter empty.  For example, 
   to NOT run a parallel compiler test:
   F77_MPI=''
   CC_MPI=''

 
4. `$ your_RunTest`
   This will execute the testing sequence.  

5. To Analyze the test results run `Analysis.py ifmpi` in the directory
   created for the compiler you are testing.  `ifmpi` is a string 
   to trigger if MPI tests are to be analyzed.  Analysis.py mpi will 
   run MPI analysis, and Analysis.py serial will run serial only 
   analysis.

The logfiles from each run will be placed in a sub-directory under 
nek5_svn/tests.  These directories are named according to the compiler
used in the test.  This parameter should be set in the your_RunTest
script.  If this is not the 
first time you have ran the Buildbot tests, you may want to move 
or delete previous logfiles to prevent false positives in the Analysis.

Analysis.py is the script that will actually run the analysis on the 
tests done in RunTests.  An ' F ' is shown for tests that have failed.

A quick way to check what tests failed is to:
   `$ Analysis.py mpi | grep ' F '`


### Overview of the Buildbot Scripts

#### RunTests:
RunTests is the main driver of the buildbot tests.  This script compiles
the nek tools, edits SIZE and .rea files, and calls ExTest and ExTestmpi

For all examples, the .map files are removed and generated from the 
Nek5000 tool, genmap. 

#### ExTestmpi:
This is the script that RunTests calls that has the set of parallel
tests for each example using the parallel compiler provided by the
F77_MPI and CC_MPI parameters.

#### ExTest:  
This is the script RunTests calls that has the set of serial tests for
each example using the compiler set by F77_SRL and CC_SRL parameters.

#### Analysis.py ifmpi	
Python script used to analyze the results of RunTests:
        -takes a string as a parameter to test mpi or not
         If the string is == "mpi" mpi analysis will be ran
         Otherwise, it will not.

	-Tests for successful 'tools' compilation
	-Tests Serial time elapsed 
	-Tests for Serial and Parallel error checks
	-Tests Examples for iteration counts in pressure solver
