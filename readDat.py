# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

f = open("gaborgen24_fMRI_Day1_103_logfile.dat", "r")
print(f.read)

CSplusShock = []
CSplusNoShock = []
GS1 = []
GS2 = []
GS3 = []

datStr = []
import pandas

with open ('gaborgen24_fMRI_Day1_103_logfile.dat') as datFile:
    data = pandas.read_csv(datFile, sep=",")
    
for line in datFile:
    print(line[12])
    datStr.append(line)
        

