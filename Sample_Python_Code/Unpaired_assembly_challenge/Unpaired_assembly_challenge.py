#!/usr/bin/env python
#Computational Genomics
#Written By Jack Zhan

#Part 1: Get and parse the reads
#Get reads from reads.fna and place them into dictionary
#All the reads come from the same synthetic genome and each is 100 nt long. 
#For simplicity, these reads don't have any quality values.


def parse_fasta(fh):
    ''' Parse FASTA into a dictionary '''
    fa = {}
    name = None
    # Part 1: compile list of lines for each sequence
    for ln in fh:
        if ln[0] == '>':  # new sequence
            name = ln[1:].split()[0]
            fa[name] = []
        else:
            # append nucleotides to current sequence
            fa[name].append(ln.rstrip())
    # Part 2: join lists into strings
    for name, nuc_list in fa.items():
        fa[name] = ''.join(nuc_list)  # join into one long string
    return fa

def make_kmer_table(seqs, k):
    ''' Given dictionary (e.g. output of parse_fasta) and integer k,
        return a dictionary that maps each k-mer to the set of names
        of reads containing the k-mer. '''
    table = {}  # maps k-mer to set of names of reads containing k-mer
    for name, seq in seqs.items():
        for i in range(0, len(seq) - k + 1):
            kmer = seq[i:i+k]
            if kmer not in table:
                table[kmer] = set()
            table[kmer].add(name)
    return table

import os
cwd = os.getcwd()

with open(cwd + "\\reads.fa") as f:
    seq_dict = parse_fasta(f)
    seq_table = make_kmer_table(seq_dict,40)

#Part 2: Build an Overlap Graph
'''
For each read AA, find the other read BB that has the longest suffix/prefix match with AA, 
i.e. a suffix of AA matches a prefix of BB. BB is AA's best buddy to the right. 
However, if there is a tie, or if the longest suffix/prefix match is less than 40 
nucleotides long, then AA has no best buddy to the right. For each read, your program 
should output either (a) nothing, if there is no best buddy to the right, or (b) a single, 
space-separated line with the IDs of AA and BB and the length of the overlap, like this:

0255/2 2065/1 88

This indicates a 88 bp suffix of the read with ID 0255/2 is a prefix of the read with ID 2065/1. 
Because of how we defined best buddy, it also means no other read besides 2065/1 has a prefix of 
88+ bp that is also a suffix of read 0255/2. A corrolary of this is that a particular read ID 
should appear in the first column of your program's output at most once. Also, since we require 
the overlap to be at least 40 bases long, no number less than 40 should every appear in the last 
column.
'''

def suffixPrefixMatch(str1, str2, min_overlap):
    ''' Returns length of longest suffix of str1 that is prefix of
        str2, as long as that suffix is at least as long as min_overlap. '''
    if len(str2) < min_overlap: return 0
    str2_prefix = str2[:min_overlap]
    str1_pos = -1
    while True:
        str1_pos = str1.find(str2_prefix, str1_pos + 1)
        if str1_pos == -1: return 0
        str1_suffix = str1[str1_pos:]
        if str2.startswith(str1_suffix): return len(str1_suffix)


output = []
score = {}
partner = {}
for name, seq in seq_dict.items():
    max_len = 0
    max_name = "NA"
    dup = False
    counter = 0
    for item in seq_table[seq[-40:]]:
        if item != name:
            match_len = suffixPrefixMatch(seq,seq_dict[item],40)
            if match_len>max_len:
                #Replace old length
                max_len = match_len
                max_name = item
                dup = False
            elif match_len == max_len:
                #Check for duplicates
                dup = True
    if dup == True:
        output.append(name+" NA 0")
        partner[name] = "NA"
        score[name+"NA"] = 0
    else: 
        output.append(name+" "+max_name+" "+str(max_len))
        partner[name] = max_name
        score[name+max_name] = max_len
    counter +=1
#Check best buddy to the left
for key, value in partner.items():
    old_key2 = ""
    for key2, value2 in partner.items():
        if value2 == key:
            if old_key2 != "":
                partner[key2] = "NA"
                score[key2+"NA"] = 0
                partner[old_key2] = "NA"
                score[old_key2+"NA"] = 0                
            old_key2 = key2

#Part 3: Build unitigs

#Build and assemble genome with best buddy algorithm 


#Function to get the previous item in the unitig
def prev(item, dictionary):
    output = "NULL"
    for key, value in dictionary.items():
        if value == item:
            output = key
    return output
    
check_list = []
output2 = []
unitigs = []
counter = 1
for item in partner:
    if item != "NA" and item not in check_list:
        forwards_item = item
        backwards_item = item
        flag1 = False
        flag2 = False
        order = []
        order.append(forwards_item)
        while True:
            #look Forewards to build unitig
            if flag1 == False:
                forwards_item = partner[forwards_item]
                if forwards_item == "NA":
                    flag1 = True
                elif forwards_item in check_list:
                    flag1 = True
                else:
                    check_list.append(forwards_item)
                    order.append(forwards_item)
            #Look Backwards to build unitig
            if flag2 == False:
                backwards_item = prev(backwards_item, partner)
                if backwards_item == "NULL":
                    flag2 = True
                else:
                    check_list.append(backwards_item)
                    order.insert(0,backwards_item)
            if flag1 == True and flag2 == True:
                break
        output2.append("START UNITIG "+str(counter) + " " + order[0])
        unitig = seq_dict[order[0]]
        for i in range(1, len(order)):
            output2.append("  "+order[i]+" "+str(score[order[i-1]+order[i]]))
            unitig += seq_dict[order[i]][score[order[i-1]+order[i]]:]
        output2.append("END UNITIG "+str(counter))
        unitigs.append(unitig)
        counter += 1

start = ""
end = ""
for unitig1 in unitigs:
    flag1 =False
    flag2 =False
    for unitig2 in unitigs:
        if unitig1 != unitig2:
            if suffixPrefixMatch(unitig1,unitig2,5)>5: flag1 = True
            if suffixPrefixMatch(unitig2,unitig1,5)>5: flag2 = True
    if flag1 == True and flag2 == False:
        #Start Unitig
        start = unitig1
    if flag2 == True and flag1 == False:
        #End Unitig
        end = unitig1
unitigs.remove(start)
unitigs.remove(end)


counter = 0
flag1 = False
used_unitig = []
added_unitig = ""
while flag1 == False:
    max_overhang = 20
    for unitig in unitigs:
        if unitig != added_unitig:
            match_len = suffixPrefixMatch(start,unitig,20)
            if match_len > max_overhang:
                max_overhang = match_len
                new_addition = unitig
    start += new_addition[match_len:]
    added_unitig = new_addition
    if new_addition not in used_unitig: used_unitig.append(new_addition)
    #All the unitigs must be used
    if len(used_unitig) == len(unitigs):
        if suffixPrefixMatch(start,end,20) >20:
            flag1 = True
            start += end[suffixPrefixMatch(start,end,20):]
    counter += 1

import sys
def write_solution(genome, per_line=60, out=sys.stdout):
    offset = 0
    out.write('>solution\n')
    while offset < len(genome):
        nchars = min(len(genome) - offset, per_line)
        line = genome[offset:offset+nchars]
        offset += nchars
        out.write(line + '\n')


write_solution(start)