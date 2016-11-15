#/usr/bin/python
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--col", "-c", type=int, default=0,  help="Spalte nach der sortiert werden soll.") 
parser.add_argument("--sep", "-s", type=int, default=0,  help="Leere Zeile alle x Zeilen einfuegen.") 
parser.add_argument("file", help="zu sortierende Datei")

args = parser.parse_args()

i = 0
for l in sorted(open(args.file,'r').readlines(), key=lambda line: float(line.strip().split()[args.col])) : 
    print l,
    i = i+1
    if args.sep > 0 :
        if i == args.sep :
            print ''
            i = 0
