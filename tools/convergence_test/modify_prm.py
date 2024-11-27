#!/usr/bin/python
# -*- coding: utf-8 -*-

# modify the parameterfile given by 'path'
# 'properties' must be a map containing the properties to change
def modify_prm(path, properties) :
    lines = open(path, 'r').readlines()
    # iterate over all lines of parameter file
    for i in range(len(lines)) :
        line = lines[i]
        # split line at '='. Before is the property
        tmp = line.split("=")
        if len(tmp) < 2:
            continue

        prop = tmp[0]
        # iterate over all properties, that must be changed and check if
        # one matches the property of the current line
        for key, value in properties.items() :
            if key == prop.strip() :
                # change property
                tmp = tmp[1].split("!")
                if len(tmp) > 1 :
                    lines[i] = "%s= %s !%s" % (prop, str(value), tmp[1])
                else :
                    lines[i] = "%s= %s\n" % (prop, str(value))
    # write parameter file
    f = open(path, 'w')
    for line in lines :
        f.write(line)
    f.close()


def read_prm(path, param) :
    lines = open(path, 'r').readlines()
    # iterate over all lines of parameter file
    for line in lines :
        # split line at '='. Before is the property
        tmp = line.split("=")
        if len(tmp) < 2:
            continue

        if tmp[0].strip() == param :
            # split tmp at '!'. Before is the  value
            tmp = tmp[1].split("!")
            if len(tmp) >= 1 :
                return tmp[0].strip()
            else :
                return None
