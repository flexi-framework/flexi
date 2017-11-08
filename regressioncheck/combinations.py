import re
import logging
import collections
import os
import tools

class Option :
    """Create an option object with a "name" and "values" similar to dict """
    def __init__(self, name, values) :
        self.name = name
        self.values = values

def splitValues(s) :
    """ split string of values at ',' but not inside brackets, since a value can be
     an array, which is written as '(/ a, b, c /)'  """
    # This is done with a regular expression. Explanation:
    #   ,       : matches comma ','
    #   \s*     : match 0 or more whitespaces (actually not necessary, since we already removed all whitespaces)
    #   (?!...) : matches if ... doesn't match next. 
    #   [^()]*  : matches all characters, except '(' or ')', 0 or more times
    #   \)      : matches closing bracket ')', the backslash is the escape-character
    return re.split(r',\s*(?![^()]*\))', s)

def isSubset(a, b) :
    """Check if the dictionary 'a' is a subset of the dictionary 'b'"""
    try :
        # build list of booleans, that contains for every key in 'a', if a[key] == b[key]
        tmp = [ a[key] == b[key] for key in a.keys() ] 
    except KeyError : # if a key of 'a' is not in 'b'
        return False 
    return all(tmp) # return True if all elements of tmp are True

def anyIsSubset(alist, b) :
    """Check if any element 'a' of the list 'alist' is a subset of the dictionary 'b'"""
    tmp = [isSubset(a, b) for a in alist] # build a list of booleans, that contains for every 'a' in alist if 'a' is a subset of 'b'
    return any(tmp)                       # return True, if any 'a' of alist is a subset of 'b'

def readKeyValueFile(filename) :
    # General worflow:
    # 1.  Read file line by line:
    # 1.1   get exclusion from line (if line starts with 'exclude:')
    # 1.2   get noCrossCombination from line (if line starts with 'nocrosscombination:')
    # 1.3   get option and its values from line ( option=value1 [,value2 [,value3 ...]] )
    found = os.path.exists(filename) # check if directory exists
    if not found :
        #raise getCombinationException(filename) # file not found
        raise Exception(tools.red("getCombination failed. file '%s' not found." % filename))

    options = []                               # list of all options
    exclusions = []                            # list of all exclusions
    noCrossCombinations = []                   # list of all noCrossCombinations

    # 1. read options and exclusions from the file
    with open(filename) as f :
        for line in f.readlines() :   # iterate over all lines of the file
            line = re.sub(r"\s+", "", line)        # remove all whitespaces 
            if line.startswith('!') : continue     # skip lines starting with a comment
            line = line.split('!')[0]              # remove comments 

            # 1.1 read an exclusion 
            if line.lower().startswith('exclude:') :
                line = line.split(':', 1)[1]       # remove everything before ':''
                ex = {}                            # new dictionary for the exclusion
                for key_value in splitValues(line):# split at ',' (but not inside brackets) and iterate over key-value-pairs 
                    (key,value) = key_value.split('=') 
                    ex[key] = value                # save key and its value in the exclusion-dictionary

                exclusions.append(ex)              # append exclusion to the list of all exclusions
                continue                           # reading of exclusion finished -> go on with next line

            # 1.2 read a noCrossCombination
            if line.lower().startswith('nocrosscombination:') :
                line = line.split(':', 1)[1]                   # remove everything before ':''
                noCrossCombination = line.split(',')           # list of keys, that should not be cross combined
                noCrossCombinations.append(noCrossCombination) # append noCrossCombination to the list of all noCrossCombinations
                continue                                       # reading of noCrossCombination finished -> go on with next line


            # 1.3 read a option and its possible values 
            if '=' in line :
                (key,values) = line.split('=',1)         # split line at '=' 
                option = Option(key,splitValues(values)) # generate new Option with a list of values (splitted at ',' but not inside brackets)
                options.append(option)                   # append option to options list, where 
                continue                                 # reading of option finished -> go on with next line

    options.sort(key=lambda option: len(option.values), reverse=True) # sort list in order to have the most varying option at the beginning

    return options, exclusions, noCrossCombinations

def getCombinations(filename) :
    # 1. get the key-value list from file
    # 1.1   get exclusion from line (if line starts with 'exclude:')
    # 1.2   get noCrossCombination from line (if line starts with 'nocrosscombination:')
    # 1.3   get option and it values from line ( option=value1 [,value2 [,value3 ...]] )
    options, exclusions, noCrossCombinations = readKeyValueFile(filename)

    # 2.  Compute combinations:
    # 2.1   count total number of all combinations
    # 2.2   build only the valid combinations (that do NOT match any exclusion)
    # 2. compute combinations
    # 2.1 count total number of possible combinations without the exclusions
    combinations = []                          # list of all VALID combinations

    noCombinationsTotal = 1
    for option in options :
        option.base = noCombinationsTotal         # save total  number of combinations of all options before this option
        noCombinationsTotal = noCombinationsTotal * len(option.values)

    logging.getLogger('logger').debug("  Total number of combinations for '%s' = %d" % (filename, noCombinationsTotal))

    if noCombinationsTotal > 10000:
        print tools.red("more than 10000 combinations in parameter.ini not allowed!")
        exit(1) # TODO: raise exception here!

    # 2.2 build all valid combinations (all that do not match any exclusion)
    for i in range(noCombinationsTotal) :         # iterate index 'i' over noCombinationsTotal
        combination = collections.OrderedDict()
        digits = collections.OrderedDict()
        # build i-th combination by adding all options with their name and a certain value
        for option in options :
            # compute index in the list of values of the option
            # Explanation with Example: 
            #   Assume you are reading the following file:
            #       opt1 = black,white
            #       opt2 = a,b,c
            #       opt3 = cat,dog,bird,snake
            #       opt4 = car,train
            #   Then yout get 2*3*4*2 = 48 combinations in total. Let us imagine the options (opt1, opt2, ...)
            #   as digits in a crazy number system, where opt1 is the digit with the lowest value and 
            #   opt4 with the highest value. Since we have different number of values for every digit, the base 
            #   of each digit is not as in a standard number system like the decimal a series of powers 10^0, 10^1, 10^2, ...
            #   In our example we get the following bases for the options:
            #       base of opt1 = 1   (first digit has base 1 in every number system)
            #       base of opt2 = 2   (since opt1 has two   values, we get a 2 here)
            #       base of opt3 = 6   (since opt2 has three values and base=2, we get 2*3=6.  This is the number of combinations of opt1 and opt2)
            #       base of opt4 = 24  (since opt3 has foure values and base=6, we get 6*4=24. This is the number of combinations of opt1, opt2 and opt3)
            #   Luckily we already stored the base while counting the number of all combinations before.
            #   We now can compute the index in the list of values of an option (which is the value of the respective digit in our crazy number system)
            #   by dividing the index i (which is a number in our crazy number system) by the base and modulo the number of values of the option.
            j = (i / option.base) % len(option.values)
            digits[option.name] = j

        for option in options : 
            combination[option.name] = option.values[digits[option.name]]

        # check if the combination is valid (does not match any exclusion)
        if anyIsSubset(exclusions, combination) : 
            continue # if any exclusion matches the combination, the combination is invalid => cycle and do not add to list of valid combinations

        skip = False
        for noCrossCombination in noCrossCombinations :
            if not all([digits[key] == digits[noCrossCombination[0]] for key in noCrossCombination]) :
                skip = True
                break
        if skip : continue
                


        # add valid combination 
        combinations.append(combination)

    logging.getLogger('logger').debug("  Number of valid combinations = %d" % len(combinations))
    return combinations, digits


def writeCombinationsToFile(combinations, path) : # write one set of parameters to a file, e.g., parameter.ini
    with open(path, 'w') as f :
        for key, value in combinations.items() :
            # for parameters with value 'crosscombinations' in the key-value pair, replace it with the value from 'crosscombinations'
            # example: N                 = crosscombinations
            #          N_Geo             = crosscombinations
            #          crosscombinations = 1,2,3,4,5
            if value == 'crosscombinations' :
                f.write("%s=%s\n" % (key,combinations.get('crosscombinations')))
            else :
                f.write("%s=%s\n" % (key, value))

#class getCombinationException(Exception) : # Exception for missing files, e.g., command_line.ini
    #def __init__(self, filename):
        #self.filename = filename
    #def __str__(self):
        #return tools.printr("getCombination failed. file '%s' not found." % (self.filename))



