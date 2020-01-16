#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 20:00:46 2019

@author: Ayca
"""
#pip install bitarray
#from bitarray import bitarray
import numpy as np
 
def fileRead(fasta):
    file = open(fasta,'r')
    sequences = []
    names = []
    for line in file:
        if line[0] =='>':
            line = line.replace(">","")
            names.append(line.replace("\n",""))
        elif line[0] != '>':
            line = line.replace("\n","")
            sequences.append(line)
    return names,sequences


#%%
import getopt

args = str(input())
args = args.split()
args = args[1:]

#bloomFilter --ref reference.fasta --query query.fasta --kmer 3 --bloomsize 10
opts,args = getopt.getopt(args, "r:q:k:b:", ['ref=', 'query=','kmer=','bloomsize='])
for i in opts:
    if i[0] == '--ref':
        ref = i[1]
    elif i[0] == '--query':
        query = i[1]
    elif i[0] == '--kmer':
        kmer = int(i[1])
    elif i[0] == '--bloomsize':
        bloomsize = int(i[1])
        
bloomFilter = np.zeros(bloomsize,dtype=bool)

ref_labels = fileRead(ref)[0]
references = fileRead(ref)[1]

query_labels = fileRead(query)[0]
query = fileRead(query)[1]

#find cannonicals for each references
cannonical = ''
cannonicals = []
cannonical_labels = ref_labels
for i in references:
    cannonical = ''
    for j in i:
        if j == 'A':
            cannonical += 'T'
        elif j == 'C':
            cannonical += 'G'
        elif j == 'G':
            cannonical += 'C'
        elif j == 'T':
            cannonical += 'A'
    cannonicals.append(cannonical)

#find cannonicals for each references 
cannonical = ''
cannonicals_query = []
cannonical_query_labels = query_labels
for i in query:
    cannonical = ''
    for j in i:
        if j == 'A':
            cannonical += 'T'
        elif j == 'C':
            cannonical += 'G'
        elif j == 'G':
            cannonical += 'C'
        elif j == 'T':
            cannonical += 'A'
    cannonicals_query.append(cannonical)   
    
#%%lists for hashing      
hash_ref = []
for i in references:
    if len(i)>= kmer:
        for j in range(0,len(i)-kmer+1):
            hash_ref.append(i[j:j+kmer])
for i in cannonicals:
    if len(i)>= kmer:
        for j in range(0,len(i)-kmer+1):
            hash_ref.append(i[j:j+kmer])
#remove duplicates
hash_ref = list(dict.fromkeys(hash_ref))

           
hash_search = []
hash_search_cannonical = []
for i in query:
    if len(i)>= kmer:
        for j in range(0,len(i)-kmer+1):
            hash_search.append(i[j:j+kmer])
#for i in cannonicals_query:
#    if len(i)>= kmer:
#        for j in range(0,len(i)-kmer+1):
#            hash_search.append(i[j:j+kmer])
#remove duplicates
hash_search = list(dict.fromkeys(hash_search))
#%%hashing functions
print("({}) Number of distinct {}-mers indexed in reference(cannonicals are included).".format(len(hash_ref),kmer))
print("({}) Number of distinct {}-mers scanned in query(cannonicals are NOT included)".format(len(hash_search),kmer))
##hash function 1
for i in hash_ref:
    res = 0
    for j in i:
        res += ord(j)
    bloomFilter[res%bloomsize] = True

#count = 0
#for i in bloomFilter:
#    if i == True:
#        count+=1
#print("count after first hash", count)

##hash function 2
"""
A:2
B:3
C:4
D:5
"""
for i in hash_ref:
    res = 0
    digit = kmer -1
    for j in i:
        if j == 'A':
            res += 2*pow(4,digit)
        elif j == 'C':
            res += 3*pow(4,digit)
        elif j == 'G':
            res += 4*pow(4,digit)
        elif j == 'T':
            res += 5*pow(4,digit)
        digit -=1
    bloomFilter[res%bloomsize] = True

#count = 0
#for i in bloomFilter:
#    if i == True:
#        count+=1
#print("count after second hash", count)

##hash function 3
"""
A:1
B:3
C:5
D:7
"""
for i in hash_ref:
    res = 0
    digit = kmer -1
    for j in i:
        if j == 'A':
            res += 1*pow(4,digit)
        elif j == 'C':
            res += 3*pow(4,digit)
        elif j == 'G':
            res += 5*pow(4,digit)
        elif j == 'T':
            res += 7*pow(4,digit)
        digit -=1
    bloomFilter[res%bloomsize] = True

count = 0
for i in bloomFilter:
    if i == True:
        count+=1
#print("Number of true values in Bloom Filter", count,"\n")
    
#%% search
found = 0
for i in range(0,len(hash_search)):
    hash1 = 0
    first_result = False
    hash2 = 0
    second_result = False
    hash3 = 0
    third_result = False
    digit = kmer-1
    flag = False
    for j in range(0,len(hash_search[i])):
        hash1 += ord((hash_search[i])[j])
        if (hash_search[i])[j] == 'A':
            hash2 += 2*pow(4,digit)
            hash3 += 1*pow(4,digit)
        elif (hash_search[i])[j] == 'C':
            hash2 += 3*pow(4,digit)
            hash3 += 3*pow(4,digit)
        elif (hash_search[i])[j] == 'G':
            hash2 += 4*pow(4,digit)
            hash3 += 5*pow(4,digit)
        elif (hash_search[i])[j] == 'T':
            hash2 += 5*pow(4,digit)
            hash3 += 7*pow(4,digit)
        digit -=1
    #bloomFilter[res%bloomsize] = True
    first_result = bloomFilter[hash1%bloomsize]
    second_result =bloomFilter[hash2%bloomsize]
    third_result =bloomFilter[hash3%bloomsize]
    if first_result and second_result and third_result:
        flag = True
        found+=1
print("({}) Number of distinct {}-mers from query found in the reference.".format(found,kmer))
        
    
        
        
    
    