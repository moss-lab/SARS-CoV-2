README for python scripts used in "A map of the SARS-CoV-2 RNA Structurome"

Below is a description of the python scripts which assisted in the bioinformatic analysis of our SARS-CoV-2 data.
The name of the script is given, followed by a brief description of what the script does, and finally a description of expected input and output files.


Python Script: ct_compare.py
This script will compare an input of two .ct files which are the same length for identical sequence input, but vary in analysis parameters (e.g. global genome models made with different pipelines, inclusion of probing constraints, or varying ScanFold scanning parameters).
This script checks each nucleotide position from a reference file and reports whether that nt is identically paired or unpaired in a target file.
From this, a percent consistency (or similarity) is calculated which is the total number of consistent positions between both files divided by the total number of positions.
The output of this analysis is a text file which reports the name of the input files used, the number of paired and unpaired nts in both files, a percent consistency between files, and then it reports a list of all positions which conflict between the 2 files.
An example output is listed here:

"Reference File:
Number of nucleotides paired: xxx
Number of nucleotides unpaired: xxx
Percent of nucleotides paired: xxx

Target File:
Number of nucleotides paired: xxx
Number of nucleotides unpaired: xxx
Percent of nucleotides paired: xxx

Comparison of reference file to target file:
Reference File= predicted_model.ct
Target File= known_model.ct
Number of consistent nts/pairings: xxx
Number of conflicting nts/pairings: xxx
Percent consistency between files: xxx

Conflicting BPs between the Reference and Target file
Position	reference(i)	reference(j)	target(i)	target(j)
	xxx			xxx				xxx				xxx			xxx		"

Command line to run the code: python ct_compare.py ref_file.ct target_file.ct output_name


Python Script: ct_sensitivity_PPV.py
This script calculates the sensitivity and PPV between two input .ct files which are the same length for identical sequence input, but vary in analysis parameters (e.g. global genome models made with different pipelines, inclusion of probing constraints, or varying ScanFold scanning parameters).
Here, we define PPV and sensitivity as:
PPV = (The # of consistent BPs between both .ct files)/(The total # of predicted BPs)
Sensitivity =  (The # of consistent BPs between both .ct files)/(The total # of known BPs)
Using these definitions, purely in silico ScanFold models were treated as the predicted models and reactivity informed models were treated as "known" structures.
When comparing two reactivity informed models, one was designated as the predicted and the other as the known and the corresponding PPV and Sensitivity were recorded in our tables accordingly.
Output from this script is a simple text file containing the name of the reference (i.e. predicted) and target (i.e. known) .ct input files and the corresponding PPV and sensitivity values, se example below:
"Reference File:	predicted_structure.ct
Target File:	known_structure.ct
Sensitivity=	xxx
PPV=	xxx"

Command line to run the code: python ct_sensitivity_PPV.py known_structure.ct predicted_structure.ct [output_name]


Python Script: zscore_conflict_analyzer.py
The input for this script is a final partners file (an output .txt file from ScanFold analysis) and the output from the ct_compare.py script detailed above.
The final partners file will be associated with one of the ct files used in the ct compare script analysis (i.e. resulting from the same ScanFold analysis).
This script is mostly used to compare in silico (purely computational) ScanFold results to reactivity informed models to see how different z-score regions agree to experimentally informed models.
This script will compare the Zavg values from a final partners file of one ScanFold structure model to the conflict list (first column) generated from the ct_compare.py script
The conflict list is any nt that had an alternative structure from two compared structural models (i.e. an unconstrained and reactivity constrained model comparison)
The script takes Zavg values and bins them from -2 to +2 at increments of 1 and then checks to see how many conflicting positions are associated with each bin
The output will be a text file containing the percent agreement for each zscore bin, an example is shown below:
"input final partner file:	Final_partners_ScanFold_output.txt
input conflict list:	ct_compare_output.txt

<=-2 nt avg. z-score similarity:	xxx
Number of similar nt <= -2 zs:	xxx
Number of conflicting nt <= -2 zs:	xxx

>=-2, <=-1 nt avg. z-score similarity:	xxx
Number of similar nt <= -1 zs:	xxx
Number of conflicting nt <= -1 zs:	xxx

>=-1, <=0 nt avg. z-score similarity:	xxx
Number of similar nt >= -1 zs:	xxx
Number of conflicting nt >= -1 zs:	xxx

>=0, <=1 nt avg. z-score similarity:	xxx
Number of similar nt >= -1 zs:	xxx
Number of conflicting nt >= -1 zs:	xxx

>=1, <=2 nt avg. z-score similarity:	xxx
Number of similar nt <= -1 zs:	xxx
Number of conflicting nt <= -1 zs:	xxx

>=2 nt avg. z-score similarity:	xxx
Number of similar nt <= -2 zs:	xxx
Number of conflicting nt <= -2 zs:	xxx
"

Command line to run the code: python  zscore_conflict_analyzer.py final_partners.txt ct_compare_output.txt output_name


The python scripts bp_compare.py and bp_sensitivity_PPV.py work the same as the analogous ct scripts, but the input is .bp files (a ScanFold output) instead of .ct files.
