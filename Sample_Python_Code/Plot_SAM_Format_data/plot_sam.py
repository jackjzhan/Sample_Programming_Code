#!/usr/bin/env python

#Jack Zhan
#This is python code to plot SAM format data

from __future__ import division
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
import re
import os
dir_path = os.getcwd()

# Set initial variables
sampleNumber = 0;
x_list = []
y_list = []
file_Path = []

#list of organisms
organism_list = ["Actinomyces radicidentis","Lactobacillus hokkaidonensis","Streptococcus salivarius",
                     "Streptococcus mutans","Treponema denticola","Streptococcus mitis","Porphyromonas gingivalis",
                     "Escherichia coli","Desulfobulbus propionicus","Granulicatella elegans","Gemella haemolysans"
                     ]
#Organisms map
map = {"NZ_CP014228.1" : "Actinomyces radicidentis",
        "NZ_AP014680.1" : "Lactobacillus hokkaidonensis",
        "NZ_AP014681.1" : "Lactobacillus hokkaidonensis",
        "NZ_AP014682.1" : "Lactobacillus hokkaidonensis",
        "NZ_CP009913.1" : "Streptococcus salivarius",
        "NC_004350.2" : "Streptococcus mutans",
        "NC_002967.9" : "Treponema denticola",
        "NC_013853.1" : "Streptococcus mitis",
        "NC_010729.1" : "Porphyromonas gingivalis",
        "NC_000913.3" : "Escherichia coli",
        "NC_014972.1" : "Desulfobulbus propionicus",
        "NZ_KI391971.1" : "Granulicatella elegans",
        "NZ_ACDZ02000015.1" : "Gemella haemolysans",
        "NZ_ACDZ02000014.1" : "Gemella haemolysans",
        "NZ_ACDZ02000013.1" : "Gemella haemolysans",
        "NZ_ACDZ02000012.1" : "Gemella haemolysans",
        "NZ_ACDZ02000011.1" : "Gemella haemolysans",
        "NZ_ACDZ02000010.1" : "Gemella haemolysans",
        "NZ_ACDZ02000009.1" : "Gemella haemolysans",
        "NZ_ACDZ02000008.1" : "Gemella haemolysans",
        "NZ_ACDZ02000007.1" : "Gemella haemolysans",
        "NZ_ACDZ02000006.1" : "Gemella haemolysans",
        "NZ_ACDZ02000005.1" : "Gemella haemolysans",
        "NZ_ACDZ02000004.1" : "Gemella haemolysans",
        "NZ_ACDZ02000003.1" : "Gemella haemolysans",
        "NZ_ACDZ02000002.1" : "Gemella haemolysans",
        "NZ_ACDZ02000001.1" : "Gemella haemolysans"
        }

#Get Number of input files. This code only works to a max of 3.
filenum = len(sys.argv) - 2
alignment_Algorithm = str(sys.argv[1])

for i in range(filenum):
    #Get file Path and set up our data storage list
    file_Path.append(sys.argv[i+2])
    x_list.append([])
    y_list.append([])
    for j in range(len(organism_list)):
        x_list[i].append([])
        y_list[i].append([])

for file_index in range(filenum):
    with open(file_Path[file_index], "rt") as finput:
        #Parsing the inputed SAM file
        #Using https://en.wikipedia.org/wiki/SAM_(file_format) as a reference
        for line in finput: 
            # Skip the header lines
            if not line.startswith('@'):
                col_index = line.rstrip().split()
                #Column 4 is POS. POS is set as 0 for an unmapped read without coordinate.
                #Using try cause bwa has messy data
                try:
                    if col_index[3] != 0 and col_index[2] != "*":
                        #Get organism
                        organism = map[col_index[2]]
                        #Get index of mapped organism
                        org_index = organism_list.index(organism)
                        #Get the cigar score in Column 5
                        cigar =  col_index[5]
                        #Seperate the cigar score. Use regex with (digit) followed by (letter)
                        match = re.findall(r'(\d+)(\w)', cigar)
                        numMatch = 0
                        numMismatch = 0
                        #Obtain Number of matches and mismatched
                        for item in match:
                            if item[1] == "M":
                                numMatch += int(item[0])
                            else:
                                numMismatch += int(item[0])
                        #Calculate values and add to list
                        if numMatch+numMismatch>0:
                            identity = numMatch/(numMatch + numMismatch)*100
                            roundX = int(col_index[3])
                            #Append data to list
                            x_list[file_index][org_index].append(roundX)
                            y_list[file_index][org_index].append(identity)
                except:
                    continue
#Color and sample list
color_list = ["red", "blue", "green"]
samp_list = ["SRS014692","SRS015055","SRS019120"]

#Plot the data
for i in range(len(organism_list)):
    #Setup plot
    plt.title('Fragment Recruitment Plot Aligned via \n' + alignment_Algorithm)
    plt.xlabel('Reference Genome: ' + organism_list[i])
    plt.ylabel('Percent Identity')
    plt.axis([0, 5000000, 80, 100])
    x = np.arange(10)
    ax = plt.subplot(111)
    # Generate the scatter plot
    for j in range(filenum):
        plt.scatter(x_list[j][i],y_list[j][i], c=color_list[j], s=1, edgecolors='none')
        ax.plot(x, j * x, c=color_list[j], label=samp_list[j])
    # Create legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.85])  
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.10), fancybox=True, shadow=True, ncol=3 )
    #Save the file
    filename = alignment_Algorithm + "_" + organism_list[i] + '.png'
    plt.savefig(dir_path + "/" + filename)
    plt.gcf().clear()
    plt.close()
