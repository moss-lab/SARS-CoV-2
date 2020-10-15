#bp_compare
"""
By Collin A. O'Leary, Moss Lab, Iowa State University
Comparison between two ScanFold .bp files when files are run on the same input fasta, but vary in analysis parameters (e.g. inclusion of probing constrants or varying scanning parameters)
Command line to run the code: python bp_compare.py ref_file.bp target_file.bp output_name
"""
import sys
import os

r_bpij = []				#A list which will hold the i and j values for the reference .bp file, the 'r' prefix is used throughout on variables relating to the reference file
t_bpij = []				#A list which will hold the i and j values for the target .bp file, the 't' prefix is used throughout on variables relating to the target file

filename1 = sys.argv[1]				#The first system argument, this needs to be the reference .bp file
filename2 = sys.argv[2]			    #The second system argument, this is the target .bp file for comparison to the reference
outputname = sys.argv[3]            #Specify an output file name

consistent_pairings = 0             #A running number of consistent basepairs between the files
conflicting_pairings = 0            #A running number of conflicting basepairs between the files

with open(filename1 , 'r') as rbp:			#Open reference file in read mode and refer to it as short name, then splits the file into lines
	r_lines = rbp.readlines()[7:]
    
with open(filename2 , 'r') as tbp:			#Open target file in read mode and refer to it as short name, then splits the file into lines
	t_lines = tbp.readlines()[7:]
   
for r_line in r_lines:                      #This loop splits the r-file lines and then iterates through it, appending the i and j nt values to a list
    data = r_line.split('\t')
    r_bpij.append((data[2],data[3]))

for t_line in t_lines:                      #This loop splits the t-file lines and then iterates through it, appending the i and j nt values to a list
    data = t_line.split('\t')
    t_bpij.append((data[2],data[3]))
    
if len(r_bpij) >= len(t_bpij):              #This determines the shorter list (they should all be the same, but sometimes files have an extra space at the end, causing syntax errors) and sets that length to t (which is used in range below) 
    t = len(t_bpij)
else:
    t = len(r_bpij)
    
inconsistencies = []                        #A list to hold inconsistent basepairs between files
for x in range(0, t):                       #A loop which iterates a number of times equal to the length of the shortest list
    if r_bpij[x] == t_bpij[x]:              #This checks whether a tuple holding i and j values for the reference file is equal to that of the target file
        consistent_pairings += 1            #If consistent, 1 is added to a rolling number which serves as a count
    else:
        conflicting_pairings += 1                           #If conflicting, 1 is added to a rolling number which serves as a count
        inconsistencies.append((x+1,r_bpij[x],t_bpij[x]))   #The inconsistent bps are added to a list along with thier position
              
r_paired = 0                    #A variable which keeps track of the number of paired nts in the reference .bp file
r_unpaired = 0                  #A variable which keeps track of the number of unpaired nts in the reference .bp file
for obj in r_bpij:              #A loop that iterates through the reference bp list and checks whether or not the nt is paired and reports to the appropriate variable
    if obj[0] == obj[1]:
        r_unpaired += 1        
    else:
        r_paired += 1
        
t_paired = 0                    #A variable which keeps track of the number of paired nts in the target .bp file
t_unpaired = 0                  #A variable which keeps track of the number of unpaired nts in the target .bp file
for obj in t_bpij:              #A loop that iterates through the target bp list and checks whether or not the nt is paired and reports to the appropriate variable
    if obj[0] == obj[1]:
        t_unpaired += 1
    else:
        t_paired += 1


pct_r_paired = (100*r_paired/(r_paired+r_unpaired))             #Calculates the percent of paired nts in the reference .bp file
pct_t_paired = (100*t_paired/(t_paired+t_unpaired))             #Calculates the percent of paired nts in the target .bp file
pct_consistent = (100*consistent_pairings/(consistent_pairings+conflicting_pairings))   #Calulates the nts which are consistently paired or unpaired between the files

with open(f"{outputname}.txt", "w") as file:      #Creates and opens an outputfile with specified output name
    first_header = 'Reference File:'              #A set of different headers to be used below
    second_header = 'Target File:'
    third_header = f'Comparison of reference file to target file:\nReference File= {filename1}\nTarget File= {filename2}'
    fourth_header = 'Conflicting     BPs between the Reference and Target file\nPosition\treference(i)\treference(j)\ttarget(i)\ttarget(j)\n'
    
    #Here headers and corresponding data are written to the output file
    file.write(first_header + '\n' + f"Number of nucleotides paired: {r_paired}" + '\n' + f"Number of nucleotides unpaired: {r_unpaired}" + '\n' + f"Percent of nucleotides paired: {pct_r_paired}%\n\n")	
    file.write(second_header + '\n' + f"Number of nucleotides paired: {t_paired}" + '\n' + f"Number of nucleotides unpaired: {t_unpaired}" + '\n' + f"Percent of nucleotides paired: {pct_t_paired}%\n\n")  
    file.write(third_header + '\n' + f"Number of consistent pairings: {consistent_pairings}" + '\n' + f"Number of conflicting pairings: {conflicting_pairings}" + '\n' + f"Percent consistency between files: {pct_consistent}%\n\n")
    file.write(fourth_header)
    for obj in inconsistencies:
        file.write(f'{obj[0]}\t{obj[1][0]}\t{obj[1][1]}\t{obj[2][0]}\t{obj[2][1]}\n'  ) #This yields a list of all inconsistent nts and thier positions
    
file.close()












