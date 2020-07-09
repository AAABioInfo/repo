#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 6/11/2020
@author: Aldo Amaya
"""
#### This set of the script will take in all your arg

#Import all of these
import os
import sys
import pandas as pd
import numpy as np
import argparse
from scipy import misc
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap



#Import your variables 

if __name__ == "__main__":
    import_args = argparse.ArgumentParser(description='import all arguments for this script')
    import_args.add_argument("-c",  "--ClinData", required=True, help="Place clinical data")
    import_args.add_argument("-dst",  "--distance", required=True, help="Directory of distance")
    import_args.add_argument("-dvt",  "--diversity", required=True, help="Directory of diversity")

    i_args = import_args.parse_args()

clinical_data= i_args.ClinData
distance= i_args.distance
diversity= i_args.diversity

####################################Demultiplex####################################
###################################################################################
###################################################################################
#########lets crete a library with the clinical data############## 

ClinL = pd.read_csv(clinical_data,sep='\t') 
# ClinLD=ClinL.set_index('Name').T.to_dict('list')
# ClinLR=ClinL.set_index('Name').T.to_dict('list')

# Using DataFrame.insert() to add a column 
ClinL["Averages"] = ""
ClinL["Std"] = ""
#Iterate through dataframe and determine average and std and append 
for index, row in ClinL.iterrows():
    tempCN=row['code_name']
    Ofile=diversity+"/"+tempCN+".diversity.txt"
    Temp_DVST_data = np.loadtxt(fname = Ofile)
    Temp_Average=np.mean(Temp_DVST_data)
    Temp_Average = np.around(Temp_Average, 3)
    Temp_STD=np.std(Temp_DVST_data)
    Temp_STD = np.around(Temp_STD, 3)
    ClinL.at[index,'Averages'] = Temp_Average
    ClinL.at[index,'Std'] = Temp_STD


#Write this to CSV file
ClinL.to_csv('clinical_data.stats.txt',index=None, sep='\t', mode='a')


#sort Dataframe and make a list of 2 top averages and lowest diversity
TempSort=ClinL.sort_values('Averages')
ID_Max2_Min = TempSort.iloc[[-1, -2,0]]
Isolate_CN=ID_Max2_Min['code_name']
max1=Isolate_CN.iloc[0]
max2=Isolate_CN.iloc[1]
Min=Isolate_CN.iloc[2]
ListPlot=[max1,max2,Min]

#make a dictionary to store x y coordinates 

# initialize dictionary 
Plot_Dict = {} 
  
# Lets plot your data
for prt in ListPlot:
    Dfile=distance+"/"+prt+".distance.txt"
    Tdata = np.loadtxt(fname =Dfile, delimiter=',')
    x = Tdata[:, 0]
    y = Tdata[:, 1]
    Temp_plotD=[x,y]
    Plot_Dict[prt]=Temp_plotD


#iterate through dictionary to print individual; figures
for key, value in Plot_Dict.items():
    tempName=key+" distance plot"
    TFile=key+".pdf"
    tempX=value[0]
    tempY=value[1]
    plt.scatter(tempX, tempY, color="black")
    plt.title(tempName)
    plt.xlabel("Coordinate B")
    plt.ylabel("Coordinate A")
    plt.savefig(TFile)




