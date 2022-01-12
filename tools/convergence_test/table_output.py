#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys

def print_header(columns) :
    for col in columns :
        l = len(col)
        sys.stdout.write(l*"═")
        sys.stdout.write("╦")
    sys.stdout.write('\n')

    for col in columns :
        sys.stdout.write(col)
        sys.stdout.write("║")
    sys.stdout.write('\n')

    for col in columns :
        l = len(col)
        sys.stdout.write(l*"═")
        sys.stdout.write("╬")
    sys.stdout.write('\n')

def print_values(columns, formats) :
    for col,form in zip(columns,formats) :
        sys.stdout.write(form % col)
        sys.stdout.write("║")
    sys.stdout.write('\n')
