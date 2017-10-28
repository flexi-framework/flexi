import numpy as np
from externalcommand import ExternalCommand
import analyze_functions
import combinations 
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
    analyze = [] # list
    options_list, _, _ = combinations.readKeyValueFile(path)
    
    options = {} # dict
    for option in options_list :
        if len(option.values) > 1 :
            options[option.name.lower()] = option.values    # set name to lower case
        else :
            options[option.name.lower()] = option.values[0] # set name to lower case

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
    #convtest_p_tolerance = float(options.get('analyze_convtest_p_tolerance',1e-2))
    convtest_p_rate      = float(options.get('analyze_convtest_p_rate',-1))
    # only do convergence test if convergence rate and tolerance >0
    if 0.0 <= convtest_p_rate <= 1.0:
        analyze.append(Analyze_Convtest_p(convtest_p_rate))

    # 2.4   h5diff (relative or absolute HDF5-file comparison of an output file with a reference file)
    h5diff_reference_file  = options.get('h5diff_reference_file',None)
    h5diff_file            = options.get('h5diff_file',None)
    h5diff_data_set        = options.get('h5diff_data_set',None)
    h5diff_tolerance_value = options.get('h5diff_tolerance_value',1e-5)
    h5diff_tolerance_type  = options.get('h5diff_tolerance_type','absolute')
    # only do h5diff test if all variables are defined
    if h5diff_reference_file and h5diff_file and h5diff_data_set :
        if h5diff_tolerance_type in ('absolute', 'delta', '--delta') :
            h5diff_tolerance_type = "--delta"
        elif h5diff_tolerance_type in ('relative', "--relative") :
            h5diff_tolerance_type = "--relative"
        else :
            raise Exception(tools.red("initialization of h5diff failed. h5diff_tolerance_type '%s' not accepted." % h5diff_tolerance_type))
        analyze.append(Analyze_h5diff(h5diff_reference_file, h5diff_file,h5diff_data_set, h5diff_tolerance_value, h5diff_tolerance_type))
    
    #exit(1)

    return analyze

#==================================================================================================
 
class Analyze() :
    total_errors = 0

class Analyze_L2(Analyze) :
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
                Analyze.total_errors+=1

    def __str__(self) :
        return "perform L2 error comparison with a pre-defined tolerance"

#==================================================================================================

class Analyze_Convtest_h(Analyze) :
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
            print mean
            
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
                    Analyze.total_errors+=1

        else :
            print "cannot perform conv test, because number of successful runs must equal the number of cells"
            print "length(runs) ",len(runs)
            print "length(conv) ",len(self.cells)
    def __str__(self) :
        return "perform L2 h-convergence test and compare the order of convergence with the polynomial degree"

#==================================================================================================

class Analyze_Convtest_p(Analyze) :
    def __init__(self, rate) :
        self.rate = rate

    def perform(self,runs) :

        # General workflow:
        # 1.  read the polynomial degree  for all runs
        # 2.  check if number of successful runs must be euqal the number of supplied cells
        # 2.2   get L2 errors of all runs and create np.array
        # 2.3   get number of variables from L2 error array
        # 2.4   determine order of convergence between two runs
        # 2.5   check if the order of convergence is always increasing with increasing polynomial degree
        # 2.6   determine success rate from increasing convergence
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
            
            # 2.6   determine success rate from increasing convergence
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
                    Analyze.total_errors+=1

                    #global_errors+=1
        else :
            print "cannot perform conv test, because number of successful runs must equal the number of polynomial degrees"
            print "length(runs) ",len(runs)
            print "length(p)    ",len(p)
    def __str__(self) :
        return "perform L2 p-convergence test and check if the order of convergence increases with smaller grid size"



#==================================================================================================

class Analyze_h5diff(Analyze,ExternalCommand) :
    def __init__(self, h5diff_reference_file, h5diff_file, h5diff_data_set, h5diff_tolerance_value, h5diff_tolerance_type) :
        self.reference_file   = h5diff_reference_file
        self.file             = h5diff_file
        self.data_set          = h5diff_data_set
        self.tolerance_value  = h5diff_tolerance_value
        self.tolerance_type   = h5diff_tolerance_type
        ExternalCommand.__init__(self)

    def perform(self,runs) :

        # General workflow:
        # 1.  iterate over all runs
        # 1.2   execute the command 'cmd' = 'h5diff -r --XXX [number] ref_file file DataArray'
        # 1.3   if the command 'cmd' returns a code != 0, set failed
        # 1.3.1   add failed info (for return a code != 0) to run
        # 1.3.2   set analyzes to fail (for return a code != 0)

        # 1.  iterate over all runs
        for run in runs :
            # 1.2   execute the command 'cmd' = 'h5diff -r [--type] [number] [ref_file] [file] [DataSetName]'
            print " "
            cmd = ["h5diff","-r",self.tolerance_type,str(self.tolerance_value),str(self.reference_file),str(self.file),str(self.data_set)]
            print "Running ["," ".join(cmd),"]",
            try :
                self.execute_cmd(cmd, run.target_directory,"h5diff") # run the code

                # 1.3   if the comman 'cmd' return a code != 0, set failed
                if self.return_code != 0 :
                    print "   tolerance_type     : "+self.tolerance_type
                    print "   tolernace_value    : "+str(self.tolerance_value)
                    print "   reference          : "+str(self.reference_file)
                    print "   file               : "+str(self.file)
                    print "   data_set           : "+str(self.data_set)

                    # 1.3.1   add failed info if return a code != 0 to run
                    if len(self.stdout) > 20 :
                        run.analyze_results.append("h5diff failed, self.return_code != 0")
                        for line in self.stdout[:10] : # print first 10 lines
                            print " "+line,
                        print " ... leaving out intermediate lines"
                        for line in self.stdout[-10:] : # print last 10 lines
                            print " "+line,
                    else :
                        print " "+str(self.stdout)
                        if len(self.stdout) == 1 :
                            run.analyze_results.append(str(self.stdout))

                    # 1.3.2   set analyzes to fail if return a code != 0
                    run.analyze_successful=False
                    Analyze.total_errors+=1

                    #global_errors+=1
            except Exception,ex :
                self.result=tools.red("h5diff failed."+str(ex)) # print result here, because it was not added in "execute_cmd"
                print " "+self.result

                # 1.3.1   add failed info if return a code != 0 to run
                run.analyze_results.append(tools.red("h5diff failed."+str(ex)))
                run.analyze_results.append(tools.red("try adding 'export PATH=/opt/hdf5/1.X/bin/:$PATH'"))

                # 1.3.2   set analyzes to fail if return a code != 0
                run.analyze_successful=False
                Analyze.total_errors+=1


    def __str__(self) :
        return "perform h5diff between two files: ["+str(self.file)+"] + reference ["+str(self.reference_file)+"]"





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
#   Examples(iExample)%ErrorStatus=3
# END IF
