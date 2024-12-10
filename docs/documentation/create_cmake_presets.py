# ---------------------------------------------------------------------------------------------
# translate FLEXI build configuration to CMake configurePreset and append to CMakePresets.json)
# usage: python create_cmake_presets.py <path/to/build/folder> --name <name of testcase>
# ---------------------------------------------------------------------------------------------
import argparse
import json
import os

# command line arguments: build directory, name of preset (name of testcase/tutorial)
parser = argparse.ArgumentParser(description='translate FLEXI build configuration to CMake preset')
parser.add_argument('build_dir')
parser.add_argument('-n','--name',type=str,help='name of preset')
args = parser.parse_args()

# read CMakePresets.json
CMakePresets_file = '../../CMakePresets.json'
with open(CMakePresets_file,'r') as file:
    data = json.load(file)

# declare empty dictionary for new preset
new_preset = {}
new_preset['name']           = args.name
new_preset['hidden']         = False
new_preset['generator']      = 'Unix Makefiles'
new_preset['cacheVariables'] = {}

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
            new_preset['cacheVariables'][key] = val

# check if passed preset name already exists in list of configurePresets
# print(new_preset)                       # print added preset
preset_list = [ cp['name'] for cp in data['configurePresets'] ]
if args.name in preset_list:  # overwrite existing entry
  print('WARNING: Desired preset name already exists -> overwriting existing configurePreset!')
  data['configurePresets'][ preset_list.index(args.name) ] = new_preset
else: # append new preset to list of configurePresets
  data['configurePresets'].append(new_preset)
json_object = json.dumps(data,indent=2)     # serialize

# dump json object to file, modifying existing CMakePresets.json
with open(CMakePresets_file,'w') as file:
    file.write(json_object)
