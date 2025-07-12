# ---------------------------------------------------------------------------------------------
# translate FLEXI build configuration to CMake configurePreset/buildPreset and append to CMakePresets.json
# usage:    python create_cmake_presets.py <path/to/build/folder> --name <name of testcase> --binaryDir <path to tutorial build directory>
# example:  python create_cmake_presets.py ../../build/ --name convtest_inviscid --binaryDir tutorials/convtest/build_inviscid
# ---------------------------------------------------------------------------------------------
import argparse
import json
import os

# command line arguments: build directory, name of preset (name of testcase/tutorial)
parser = argparse.ArgumentParser(description='translate FLEXI build configuration to CMake preset')
parser.add_argument('build_dir')
parser.add_argument('-n','--name',type=str,help='name of preset')
parser.add_argument('-b','--binaryDir',type=str,help='path to tutorial build directory')
args = parser.parse_args()

# read CMakePresets.json
CMakePresets_file = '../../CMakePresets.json'
with open(CMakePresets_file,'r') as file:
    data = json.load(file)

# declare empty dictionary for new presets
configure_preset = {}
configure_preset['name']           = args.name
configure_preset['binaryDir']      = args.binaryDir
configure_preset['hidden']         = False
configure_preset['generator']      = 'Unix Makefiles'
configure_preset['cacheVariables'] = {}
build_preset                    = {}
build_preset['name']            = args.name
build_preset['configurePreset'] = args.name

# read build options from CMakeCache.txt to define preset
CMakeCache_file = 'CMakeCache.txt'
with open(os.path.join(args.build_dir,CMakeCache_file),'r') as file:
    # process file line-wise
    for line in file:
        # extract configured build options
        # given in the form FLEXI_<BUILD_OPTION>:<TYPE>=<DEFAULT>, for example: FLEXI_2D:BOOL=OFF
        # note that line string ends with EOL, hence omit last character when checking with 'endswith' method
        if line.startswith('FLEXI_') or line.startswith('LIBS_BUILD_') or line.startswith('CMAKE_BUILD_TYPE') or (line.startswith('POSTI') and line[:-1].endswith('ON')):
            # skip unconfigured options (internal cache entries)
            if 'INTERNAL' in line:
                continue
            # add entry to dictionary: build option -> value
            key = line.split(':')[0]        # build option
            val = line.split('=')[-1][:-1]  # configured value (remove trailing EOL-character by extracting substring [:-1])
            configure_preset['cacheVariables'][key] = val

# check if passed preset name already exists in list of configurePresets
# print(configure_preset)                       # print added preset
preset_list = [ cp['name'] for cp in data['configurePresets'] ]
if args.name in preset_list:  # overwrite existing entry
  print('WARNING: Desired preset name already exists -> overwriting existing configurePreset and buildPreset!')
  data['configurePresets'][ preset_list.index(args.name) ] = configure_preset
  data['buildPresets'][ preset_list.index(args.name) ]     = build_preset
else: # append new preset to list of configurePresets
  data['configurePresets'].append(configure_preset)
  data['buildPresets'].append(build_preset)
json_object = json.dumps(data,indent=2)     # serialize

# dump json object to file, modifying existing CMakePresets.json
with open(CMakePresets_file,'w') as file:
    file.write(json_object)
