#!/usr/bin/env python

from csv import reader
from csv import writer
import sys, getopt

def main(argv):
   inputfile = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
   except getopt.GetoptError:
      print 'addZeros.py -i <inputfile> -o <outputfile>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'addZeros.py -i <inputfile> -o <outputfile>'
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
   print 'Input file is "', inputfile
   print 'Output file is "', outputfile

   default_text = '0'
   # Open the input_file in read mode and output_file in write mode
   with open(inputfile, 'r') as read_obj, \
           open(outputfile, 'w') as write_obj:
       # Create a csv.reader object from the input file object
       csv_reader = reader(read_obj)
       # Create a csv.writer object from the output file object
       csv_writer = writer(write_obj)
       # Read each row of the input csv file as list
       for row in csv_reader:
           # Append the default text in the row / list
           row.append(default_text)
           # Add the updated row / list to the output file
           csv_writer.writerow(row)

if __name__ == "__main__":
   main(sys.argv[1:])







