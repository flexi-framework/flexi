import re

class Option :
    def __init__(self, name, values) :
        self.name = name
        self.values = values

def splitValues(s) :
    # split string of values at ',' but not inside brackets, since a value can be
    # an array, which is written as '(/ a, b, c /)' 
    # This is done with a regular expression. Explanation:
    #   ,       : matches comma ','
    #   \s*     : match 0 or more whitespaces (actually not necessary, since we already removed all whitespaces)
    #   (?!...) : matches if ... doesn't match next. 
    #   [^()]*  : matches all characters, except '(' or ')', 0 or more times
    #   \)      : matches closing bracket ')', the backslash is the escape-character
    return re.split(r',\s*(?![^()]*\))', s)

def getCombinations(filename) :
    # General worflow:
    # 1.  Read file line by line:
    # 1.1   get exclusion from line (if line starts with 'exclude:')
    # 1.2   get option and it values from line ( option=value1 [,value2 [,value3 ...]] )
    # 2.  Compute combinations:
    # 2.1   count total number of all combinations
    # 2.2   build only the valid combinations (that do NOT match any exclusion)

    options = []                               # list of all options
    exclusions = []                            # list of all exclusions
    combinations = []                          # list of all VALID combinations

    # 1. read options and exclusions from the file
    for line in open(filename).readlines() :   # iterate over all lines of the file
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

        # 1.2 read a option and its possible values 
        if '=' in line :
            (key,values) = line.split('=',1)         # split line at '=' 
            option = Option(key,splitValues(values)) # generate new Option with a list of values (splitted at ',' but not inside brackets)
            options.append(option)                   # append option to options list, where 
            continue                                 # reading of option finished -> go on with next line

    # 2. compute combinations
    # 2.1 count total number of possible combinations without the exclusions
    noCombinationsTotal = 1
    for option in options :
        option.base = noCombinationsTotal         # save total  number of combinations of all options before this option
        noCombinationsTotal = noCombinationsTotal * len(option.values)

    print "  Total number of combinations for '%s' = " % filename, noCombinationsTotal

    # 2.2 build all valid combinations (all that do not match any exclusion)
    for i in range(noCombinationsTotal) :         # iterate index 'i' over noCombinationsTotal
        combination = {}
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
            combination[option.name] = option.values[j]

        # check if valid the combination is valid (does not match any exclusion)
        doExclude = False
        for ex in exclusions :     # iterate over all exclusions and check if any exclusion matches the actual combination
            allMatch = True        # assume that the values of all keys of the exclusion match with the combination
            for name in ex.keys() :# iterate over all names of the exclusion and check if all their values match with the respective value in the combination
                allMatch = allMatch and combination[name] == ex[name]
            if allMatch :          # all values of the exclusion match with the combination => exclusion matches => INVALID
                doExclude = True
                break
                
        if doExclude : continue # combination matches exclude => cycle and do not add to list of valid combinations

        # add valid combination 
        combinations.append(combination)

    print "  Number of valid combinations =", len(combinations)
    return combinations

