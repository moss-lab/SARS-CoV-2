#bp_sensitivity_PPV 
"""
By Collin A. O'Leary, Moss Lab, Iowa State University
Sensitivity and PPV comparison between two ScanFold .bp files when files are run on the same input fasta, but vary in analysis parameters (e.g. inclusion of probing constrants or known structures)
Command line to run the code: python bp_sensitivity_PPV.py ref_file.bp target_file.bp output_name
"""
import sys
import os

r_bpij = []				#A list which will hold the paired i and j values for the reference .bp file, the 'r' prefix is used throughout on variables relating to the reference file
t_bpij = []				#A list which will hold the paired i and j values for the target .bp file, the 't' prefix is used throughout on variables relating to the target file

filename1 = sys.argv[1]				#The first system argument, this needs to be the reference file (in silico)
filename2 = sys.argv[2]			    #The second system argument, this is the target file for comparison to the reference (constrained)
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
    
if len(r_bpij) >= len(t_bpij):              #This determines the shorter list and sets that length to t (which is used in range below) 
    t = len(t_bpij)
else:
    t = len(r_bpij)
    
for x in range(0, t):                                                   #A loop which iterates a number of times equal to the length of the shortest list
    if r_bpij[x][0] != r_bpij[x][1] and r_bpij[x] == t_bpij[x]:         #This checks whether a tuple holding i and j values for the reference file is paired equal to that of the target file
        consistent_pairings += 1                                        #If consistent, 1 is added to a rolling number which serves as a count
    elif r_bpij[x][0] != r_bpij[x][1] and r_bpij[x] != t_bpij[x]:
        conflicting_pairings += 1                                       #If conflicting, 1 is added to a rolling number which serves as a count
    else:                                                               #Unpaired nts are skipped
        pass
        
r_paired = 0                    #A variable which keeps track of the number of paired nts in the reference .bp file
r_unpaired = 0                  #A variable which keeps track of the number of unpaired nts in the target .bp file
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

sensitivity = consistent_pairings/t_paired      #calculates the sensitivity

PPV = consistent_pairings/r_paired              #Calculates PPV


with open(f"{outputname}.txt", "w") as file:      #Creates and opens an outputfile with specified output name
    first_header = 'Reference (Predicted) File:'              #A set of different headers to be used below
    second_header = 'Target (Known) File:'
    
    #Here headers and corresponding data are written to the output file
    file.write(f'{first_header} {filename1}\n{second_header} {filename2}\nSensitivity={sensitivity}\nPPV={PPV}')
file.close()
print(sensitivity)
print(PPV)              #Sensitivity and PPV are printed as well as being written to an output file












