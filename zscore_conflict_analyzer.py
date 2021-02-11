#zscore_conflict_analyzer
'''
By Collin A. O'Leary, Moss Lab, Iowa State University
This script will compare the Zavg values from a final partners file of one ScanFold structure model to the conflict list generated from the ct_compare.py script
The conflict list is any nt that had an alternitive structure from two compared structural models (i.e. an unconstrained and reactivity constrained model comparison)
Zavg values are binned from -2 to +2 at increments of 1 and then further divided based on their conflicts
The output will be a text file containing the percent agreement for each zscore bin
To use: python  zscore_conflict_analyzer.py final_partners.txt ct_compare_output.txt output_name
'''

import sys
import os

filename1 = sys.argv[1]				#This needs to be final partners file from ScanFold output
filename2 = sys.argv[2]			    #This is the output file from the ct_compare.py script... verify the slicing starts at the right index below
output_name = sys.argv[3]           #Need to specify an output name here, no spaces allowed



with open(filename1 , 'r') as fp:			#Open final partners file and seperate it by lines
	fp_lines = fp.readlines()[1:]  

with open(filename2 , 'r') as cp:			#Open ct_compare.py output file and seperate it by lines starting at the confliclt list
	cp_lines = cp.readlines()[19:] 

conflicts = []                      #A list to hold the conflicting nt position
for cp_line in cp_lines:            #A loop to write the conflicting position to the list
    data = cp_line.split()
    conflicts.append(int(data[0]))

ntzs = []                           #A list to hold the nt position and associate Zavg value from the final partners file
x = 1                               #Used to keep track of nt position
for fp_line in fp_lines:            #This loop writes the positin and Zavg values to the list
    data = fp_line.split()
    ntzs.append((x,float(data[4])))
    x += 1

zsn2_sim = []                       #A list to hold the position of nts with Zavg <-2 that are in agreement
zsn2_con = []                       #A list to hold the position of nts with Zavg <-2 that are in conflict
zsn1_sim = []                       #A list to hold the position of nts with Zavg >-2, <-1 that are in agreement
zsn1_con = []                       #A list to hold the position of nts with Zavg >-2, <-1 that are in conclict
zsn0_sim = []                       #A list to hold the position of nts with Zavg >-1, <0 that are in agreement
zsn0_con = []                       #A list to hold the position of nts with Zavg >-1, <0 that are in conclict

zsp2_sim = []                       #A list to hold the position of nts with Zavg >2 that are in agreement
zsp2_con = []                       #A list to hold the position of nts with Zavg >2 that are in conflict
zsp1_sim = []                       #A list to hold the position of nts with Zavg <2, >1 that are in agreement
zsp1_con = []                       #A list to hold the position of nts with Zavg <2, >1 that are in conflict
zsp0_sim = []                       #A list to hold the position of nts with Zavg <1, >0 that are in agreement
zsp0_con = []                       #A list to hold the position of nts with Zavg <1, >0 that are in conflict

for x in ntzs:                  #A loop that iterates through the ntzs list, and first determines what the Zavg value is, then checks whether its corresponding position is in conflict. It then adds the position to the appropriate bin list
    if x[1] <= -2:
        if x[0] in conflicts:
            zsn2_con.append(x[0])
        else:
            zsn2_sim.append(x[0])
    elif x[1] <= -1:
        if x[0] in conflicts:
            zsn1_con.append(x[0])
        else:
            zsn1_sim.append(x[0])
    elif x[1] <= 0:
        if x[0] in conflicts:
            zsn0_con.append(x[0])
        else:
            zsn0_sim.append(x[0])
    elif x[1] >= 2:
        if x[0] in conflicts:
            zsp2_con.append(x[0])
        else:
            zsp2_sim.append(x[0])
    elif x[1] >= 1:
        if x[0] in conflicts:
            zsp1_con.append(x[0])
        else:
            zsp1_sim.append(x[0])
    elif x[1] >= 0:
        if x[0] in conflicts:
            zsp0_con.append(x[0])
        else:
            zsp0_sim.append(x[0])
    else:
        print('error!!')
        
if (len(zsn2_sim)+len(zsn2_con)) == 0:          #This checks whether there are any values in either bin list, if there are no values then N/A, else it calculates the percent agreement of that bin
    n2percent = "N/A"
else:
    n2percent = ((len(zsn2_sim)/(len(zsn2_sim)+len(zsn2_con)))*100) 
    
if (len(zsn1_sim)+len(zsn1_con)) == 0:          #Same as above
    n2n1percent = "N/A"
else:
    n2n1percent = ((len(zsn1_sim)/(len(zsn1_sim)+len(zsn1_con)))*100)   

if (len(zsn0_sim)+len(zsn0_con)) == 0:          #Same as above
    n1n0percent = "N/A"
else:
    n1n0percent = ((len(zsn0_sim)/(len(zsn0_sim)+len(zsn0_con)))*100)
    
if (len(zsp0_sim)+len(zsp0_con)) == 0:          #Same as above
    p0p1percent = "N/A"
else:
    p0p1percent = ((len(zsp0_sim)/(len(zsp0_sim)+len(zsp0_con)))*100)
    
if (len(zsp1_sim)+len(zsp1_con)) == 0:          #Same as above
    p1p2percent = "N/A"
else:
    p1p2percent = ((len(zsp1_sim)/(len(zsp1_sim)+len(zsp1_con)))*100)
    
if (len(zsp2_sim)+len(zsp2_con)) == 0:          #Same as above
    p2percent = "N/A"
else:
    p2percent = ((len(zsp2_sim)/(len(zsp2_sim)+len(zsp2_con)))*100)

        
with open(f"{output_name}.txt", "w") as file:      #Creates and opens an outputfile with the specified output name, then writes the calculated metrics along with input file names
    file.write(f'input final partner file:\t{filename1}\n')
    file.write(f'input conflict list:\t{filename2}\n\n')
    file.write(f'<=-2 nt avg. z-score similarity:\t{n2percent}%\n# of similar nt <= -2 zs:\t{len(zsn2_sim)}\n# of conflicting nt <= -2 zs:\t{len(zsn2_con)}\n\n')
    file.write(f'>=-2, <=-1 nt avg. z-score similarity:\t{n2n1percent}%\n# of similar nt >=-2 and <=-1 zs:\t{len(zsn1_sim)}\n# of conflicting nt >=-2 and <=-1 zs:\t{len(zsn1_con)}\n\n')
    file.write(f'>=-1, <=0 nt avg. z-score similarity:\t{n1n0percent}%\n# of similar nt >=-1 and <=0 zs:\t{len(zsn0_sim)}\n# of conflicting nt >=-1 and <=0 zs:\t{len(zsn0_con)}\n\n')
    file.write(f'>=0, <=1 nt avg. z-score similarity:\t{p0p1percent}%\n# of similar nt >=0 and <=1 zs:\t{len(zsp0_sim)}\n# of conflicting nt >=0 and <=1 zs:\t{len(zsp0_con)}\n\n')
    file.write(f'>=1, <=2 nt avg. z-score similarity:\t{p1p2percent}%\n# of similar nt >=1 and <=2 zs:\t{len(zsp1_sim)}\n# of conflicting nt >=1 and <=2 zs:\t{len(zsp1_con)}\n\n')
    file.write(f'>=2 nt avg. z-score similarity:\t{p2percent}%\n# of similar nt >= 2 zs:\t{len(zsp2_sim)}\n# of conflicting nt >= 2 zs:\t{len(zsp2_con)}\n\n')
file.close()