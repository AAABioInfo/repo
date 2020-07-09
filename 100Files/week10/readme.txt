# This readme will tell you how to run the pipeline that will do the following steps. 

This was run on the following:
Python 3.6.10 
Numpy 1.18.4

In order to run you will need the following: 
-clinical data (tagged with -c)
-directory with distance files (tagged with -dst)
-directory with diversity Scores (tagged with -dvt)

To run this script with our specific data you will execute as below: 
python3 pipeline.py -c clinical_data.txt -dst distanceFiles -dvt diversityScores




Script performs the tast below:

1) Read in the data clinical file (clinical_data.txt) as a dataframe 
and create two new columns called averages and std. 
--Done--

You will read in each sample's diversity score list (found in diversityScores) 
and append the mean value and the standard deviation of the sample's diversity score. 

--Done--

Then output the new clinical data (which includes the original clinical information 
from clinical_data.txt along with the two new columns to a file 
called clinical_data.stats.txt.)

--Done--

This should be done using pandas and numpy python package.
--Done--


2) Then find the animals (code names) that correspond with the two highest average 
diversity scores and one with the lowest average diversity score (three animals in 
total).  

--Done--


Find their respective txt file inside the distanceFiles directory and plot 
and save a scatter plot for each animal (three plots in total, and can be saved as 
png files). 

--Done--

Add a title for each plot that has the animal's name in it. 

--Done--


