import numpy as np
from loop import Loop
import analyze_functions
import combinations 
import check
import tools

#==================================================================================================

def getAnalyzes(path, example) :

    # General workflow:
    # 1.  Read the analyze options from file 'path' into dict 'options'
    # 2.  Initialize analyze functions
    # 2.1   L2 error
    # 2.2   h-convergence test
    # 2.3   p-convergence test
    # 2.4   h5diff (relative or absolute HDF5-file comparison of an output file with a reference file)

    # 1.  Read the analyze options from file 'path'
    analyze = []
    options_list, _, _ = combinations.readKeyValueFile(path)
    
    options = {}
    for option in options_list :
        if len(option.values) > 1 :
            options[option.name.lower()] = option.values
        else :
            options[option.name.lower()] = option.values[0]

    # 2.1   L2 error
    L2_tolerance = float(options.get('analyze_l2',-1.))
    if L2_tolerance > 0 :
        analyze.append(Analyze_L2(L2_tolerance))
    
    # 2.2   h-convergence test
    convtest_h_cells     = [float(cell) for cell in options.get('analyze_convtest_h_cells',['-1.'])]
    convtest_h_tolerance = float(options.get('analyze_convtest_h_tolerance',1e-2))
    convtest_h_rate      = float(options.get('analyze_convtest_h_rate',1))
    # only do convergence test if supplied cells count > 0
    if min(convtest_h_cells) > 0 and convtest_h_tolerance > 0 and 0.0 <= convtest_h_rate <= 1.0:
        analyze.append(Analyze_Convtest_h(convtest_h_cells,convtest_h_tolerance,convtest_h_rate))

    # 2.3   p-convergence test
    convtest_p_tolerance = float(options.get('analyze_convtest_p_tolerance',1e-2))
    convtest_p_rate      = float(options.get('analyze_convtest_p_rate',-1))
    # only do convergence test if convergence rate and tolerance >0
    if convtest_p_tolerance > 0 and 0.0 <= convtest_p_rate <= 1.0:
        analyze.append(Analyze_Convtest_p(convtest_p_tolerance,convtest_p_rate))

    # 2.4   h5diff (relative or absolute HDF5-file comparison of an output file with a reference file)
    h5diff_reference_file = options.get('h5diff_reference_file',None)
    h5diff_file = options.get('h5diff_file',None)
    h5diff_name = options.get('h5diff_name',None)
    # only do h5diff test if all variables are defined
    if h5diff_reference_file and h5diff_file and h5diff_name :
        print h5diff_reference_file
        print h5diff_file
        print h5diff_name
        analyze.append(Analyze_h5diff(h5diff_reference_file,h5diff_file,h5diff_name))
    
    #exit(1)

    return analyze

#==================================================================================================
 
class Analyze_L2() :
    def __init__(self, L2_tolerance) :
        self.L2_tolerance = L2_tolerance

    def perform(self,runs) :

        # General workflow:
        # 1.  Iterate over all runs
        # 1.1   read L2 errors from 'std.out' file
        # 1.2   if one L2 errors is larger than the tolerance -> fail
        # 1.3   append info for summary of errors
        # 1.4   set analyzes to fail

        # 1.  Iterate over all runs
        for run in runs :
            
            # 1.1   read L2 errors from 'std.out' file
            L2_errors = np.array(analyze_functions.get_last_L2_error(run.stdout))
            
            # 1.2   if one L2 errors is larger than the tolerance -> fail
            if (L2_errors > self.L2_tolerance).any() :
                print tools.red("analysis failed: L2 error >"+str(self.L2_tolerance))
                
                # 1.3   append info for summary of errors
                run.analyze_results.append("analysis failed: L2 error >"+str(self.L2_tolerance))
                #global_errors+=1

                # 1.4   set analyzes to fail
                run.analyze_successful=False
                run.total_errors+=1

    def __str__(self) :
        return "perform L2 error comparison with a pre-defined tolerance"

#==================================================================================================

class Analyze_Convtest_h() :
    def __init__(self, cells, tolerance, rate) :
        self.cells = cells
        self.tolerance = tolerance
        self.rate = rate

    def perform(self,runs) :

        # General workflow:
        # 1.  check if number of successful runs is euqal the number of supplied cells
        # 1.1   read the polynomial degree from the first run -> must not change!
        # 1.2   get L2 errors of all runs and create np.array
        # 1.3   get number of variables from L2 error array
        # 1.4   determine order of convergence between two runs
        # 1.5   determine success rate by comparing the relative convergence error with a tolerance
        # 1.6   compare success rate with pre-defined rate
        # 1.7     interate over all runs
        # 1.7.1   add failed info if success rate is not reached to all runs
        # 1.7.2   set analyzes to fail if success rate is not reached for all runs

        # 1.  check if number of successful runs is euqal the number of supplied cells
        if len(self.cells) == len(runs) :

            # 1.1   read the polynomial degree from the first run -> must not change!
            p = float(runs[0].parameters.get('N',-1))

            # 1.2   get L2 errors of all runs and create np.array
            all_L2_errors = np.array([analyze_functions.get_last_L2_error(run.stdout) for \
                    run in runs])
            all_L2_errors = np.transpose(all_L2_errors)

            # 1.3   get number of variables from L2 error array
            nVar = len(all_L2_errors)
            print tools.blue("L2 errors nVar="+str(nVar))
            print all_L2_errors

            # 1.4   determine order of convergence between two runs
            all_L2_order = \
                    np.array([analyze_functions.calcOrder_h(self.cells,all_L2_errors[x]) for \
                    x in range(nVar)])
            print tools.blue("L2 orders for nVar="+str(nVar))
            print all_L2_order

            # determine average convergence rate
            mean = [np.mean(all_L2_order[x]) for x in range(nVar)]
            
            # 1.5   determine success rate by comparing the relative convergence error with a tolerance
            print "relative order error ",[abs(mean[x]/(3+1)-1) for x in range(nVar)]
            success = [abs(mean[x]/(p+1)-1) < self.tolerance for x in range(nVar)]
            print "success convergence: ",success

            # 1.6   compare success rate with pre-defined rate, fails if not reached
            if float(sum(success))/nVar >= self.rate :
                print tools.blue("h-convergence successful")
            else :
                print tools.red("h-convergence failed"+"\n"+\
                        "success rate="+str(float(sum(success))/nVar)+\
                        " tolerance rate="+str(self.rate))

                # 1.7     interate over all runs
                for run in runs :

                    # 1.6.1   add failed info if success rate is not reached to all runs
                    run.analyze_results.append("analysis failed: h-convergence "\
                            +str(success))

                    # 1.6.2   set analyzes to fail if success rate is not reached for all runs
                    run.analyze_successful=False
                    run.total_errors+=1

        else :
            print "cannot perform conv test, because number of successful runs must equal the number of cells"
            print "length(runs) ",len(runs)
            print "length(conv) ",len(self.cells)
    def __str__(self) :
        return "perform L2 h-convergence test and compare the order of convergence with the polynomial degree"

#==================================================================================================

class Analyze_Convtest_p() :
    def __init__(self, tolerance, rate) :
        self.tolerance = tolerance
        self.rate = rate

    def perform(self,runs) :

        # General workflow:
        # 1.  read the polynomial degree  for all runs
        # 2.  check if number of successful runs must be euqal the number of supplied cells
        # 2.2   get L2 errors of all runs and create np.array
        # 2.3   get number of variables from L2 error array
        # 2.4   determine order of convergence between two runs
        # 2.5   check if the order of convergence is always increasing with increasing polynomial degree
        # 2.6   determine success rate by comparing the relative convergence error with a tolerance
        # 2.7   compare success rate with pre-defined rate, fails if not reached 
        # 2.8   interate over all runs
        # 2.8.1   add failed info if success rate is not reached to all runs
        # 2.8.1   set analyzes to fail if success rate is not reached for all runs

        # 1.  read the polynomial degree  for all runs
        p = [float(run.parameters.get('N',-1)) for run in runs] # get polynomial degree

        # 2   check if number of successful runs must be euqal the number of supplied cells
        if len(p) == len(runs) :

            # 2.2   get L2 errors of all runs and create np.array
            all_L2_errors = np.array([analyze_functions.get_last_L2_error(run.stdout) for \
                    run in runs])
            all_L2_errors = np.transpose(all_L2_errors)

            # 2.3   get number of variables from L2 error array
            nVar = len(all_L2_errors)
            
            print tools.blue("L2 errors nVar="+str(nVar))
            print all_L2_errors

            # 2.4   determine order of convergence between two runs
            all_L2_order = \
                    np.array([analyze_functions.calcOrder_p(p,all_L2_errors[x]) for \
                    x in range(nVar)])
            print tools.blue("L2 orders for nVar="+str(nVar))
            print all_L2_order

            # 2.5   check if the order of convergence is always increasing with increasing polynomial degree
            increasing = []
            for j in range(nVar) :
                increasing_run = []
                for i in range(1,len(p)) :
                    increasing_run.append(p[i]>p[i-1])
                #print increasing_run
                increasing.append(all(increasing_run))
            print tools.blue("Increasing order of convergence")
            print tools.blue(str(increasing))
            
            # 2.6   determine success rate by comparing the relative convergence error with a tolerance
            success = [increasing[x] for x in range(nVar)]
            print "success convergence: ",success

            # 2.7   compare success rate with pre-defined rate, fails if not reached
            if float(sum(success))/nVar >= self.rate :
                print tools.blue("p-convergence successful")
            else :
                print tools.red("p-convergence failed"+"\n"+\
                        "success rate="+str(float(sum(success))/nVar)+\
                        " tolerance rate="+str(self.rate))

                # 2.8     interate over all runs
                for run in runs :

                    # 2.8.1   add failed info if success rate is not reached to all runs
                    run.analyze_results.append("analysis failed: p-convergence "+str(success))

                    # 2.8.2   set analyzes to fail if success rate is not reached for all runs
                    run.analyze_successful=False
                    run.total_errors+=1

                    #global_errors+=1
        else :
            print "cannot perform conv test, because number of successful runs must equal the number of polynomial degrees"
            print "length(runs) ",len(runs)
            print "length(p)    ",len(p)
    def __str__(self) :
        return "perform L2 p-convergence test and check if the order of convergence increases with smaller grid size"



#==================================================================================================

class Analyze_h5diff(Loop) :
    def __init__(self, h5diff_reference_file, h5diff_file, h5diff_name) :
        self.reference_file = h5diff_reference_file
        self.file           = h5diff_file
        self.name           = h5diff_name
        Loop.__init__(self, None, 'h5diff', -1, mkdir=False)

    def perform(self,runs) :
        print self.reference_file
        print self.file
        print self.name

        # General workflow:
        # 1.  iterate over all runs
        # 1.2   select relative or absolute comparison
        # 1.2   set directory in which the program is executed
        # 1.3   execute the command 'cmd' = 'h5diff -r --XXX [number] ref_file file DataArray'
        # 1.4   if the comman 'cmd' return a code != 0, set failed
        # 1.4.1   add failed info if return a code != 0 to run
        # 1.4.2   set analyzes to fail if return a code != 0

        # 1.  iterate over all runs
        for run in runs :
            h5diff = "/opt/hdf5/1.8.18/bin/h5diff"

            # 1.2   select relative or absolute comparison
            if 1 == 1 :
                cmd = [h5diff,"-r","--delta","1e-5",str(self.reference_file),str(self.file),str(self.name)]
            else :
                cmd = [h5diff,"-r","--relative","1e-5",str(self.reference_file),str(self.file),str(self.name)]
            print "Running ["," ".join(cmd),"]",

            # 1.2   set directory in which the program is executed
            self.target_directory = run.target_directory

            # 1.3   execute the command 'cmd' = 'h5diff -r --XXX [number] ref_file file DataArray'
            self.execute_cmd(cmd) # run the code

            # 1.4   if the comman 'cmd' return a code != 0, set failed
            if self.return_code != 0 :

                # 1.4.1   add failed info if return a code != 0 to run
                run.analyze_results.append("h5diff failed")

                # 1.4.2   set analyzes to fail if return a code != 0
                run.analyze_successful=False
                run.total_errors+=1

                #global_errors+=1

    def __str__(self) :
        return "perform h5diff between two files"





# ! set ErrorStatus
# IF(iSTATUS.EQ.0)THEN
#   RETURN ! all is safe
# ELSEIF(iSTATUS.EQ.-5)THEN
#   SWRITE(UNIT_stdOut,'(A)')  ' h5diff: arrays in h5-files have different ranks.'
#   Examples(iExample)%ErrorStatus=5
# ELSEIF(iSTATUS.EQ.2)THEN
#   SWRITE(UNIT_stdOut,'(A)')  ' h5diff: file to compare not found.'
#   Examples(iExample)%ErrorStatus=5
# ELSEIF(iSTATUS.EQ.127)THEN
#   SWRITE(UNIT_stdOut,'(A)')  ' h5diff executable could not be found.'
#   Examples(iExample)%ErrorStatus=5
# ELSE!IF(iSTATUS.NE.0) THEN
#   SWRITE(UNIT_stdOut,'(A)')  ' HDF5 Datasets do not match! Error in computation!'
#   SWRITE(UNIT_stdOut,'(A)')  '    Type               : '//ADJUSTL(TRIM(Examples(iExample)%H5diffToleranceType))
#   SWRITE(UNIT_stdOut,'(A)')  '    tmpTol             : '//ADJUSTL(TRIM(tmpTol))
#   SWRITE(UNIT_stdOut,'(A)')  '    H5DIFF             : '//ADJUSTL(TRIM(H5DIFF))
#   SWRITE(UNIT_stdOut,'(A)')  '    ReferenceStateFile : '//TRIM(Examples(iExample)%H5DIFFReferenceStateFile)
#   SWRITE(UNIT_stdOut,'(A)')  '    CheckedFileName    : '//TRIM(Examples(iExample)%H5DIFFCheckedStateFile)
#   Examples(iExample)%ErrorStatus=3
# END IF
