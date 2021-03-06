#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 15:32:31 2020

@author: zhanf
"""

import numpy as np
import math
import mmh3
from bitarray import bitarray

class BloomFilter():
   
    def __init__(self, k_mers, fp_rate):
        self.num_of_kmers = len(k_mers)
        self.bloomsize = self.bloomsizecalc(len(k_mers), fp_rate)
        self.hash_count = self.hashfunctions(len(k_mers), self.bloomsize)    
        self.bits = bitarray(self.bloomsize)
        self.bits.setall(0)
        
        
        
    def bloomsizecalc(self, num_of_kmers, fp_rate):
        denom = (math.log(2))**2
        numer = num_of_kmers*math.log(fp_rate)
        return int(-(numer/denom))
    
    def hashfunctions(self, num_of_kmers, bloom):
        return int((bloom/num_of_kmers)*math.log(2))
    
    def add(self, k_mer):
        for i in range(self.hash_count):
            n = mmh3.hash(k_mer, i) % self.bloomsize
            self.bits[n] = 1
            
    def check(self, k_mer):
        for i in range(self.hash_count):
            n = n = mmh3.hash(k_mer, i) % self.bloomsize
            if self.bits[n] == False:
                return False
        return True
        
class chromosome():
    
    def __init__(self, seqobject, k_length):
        self.name = seqobject.id
        self.sequence = str(seqobject.seq)
        self.kmers = self.get_kmers(self.sequence, k_length)
        self.length = len(self.sequence)
    
    def get_kmers(self, sequence, k_length):
        k_mers = []
        for i in range(len(sequence)-k_length):
            k_mers.append(sequence[i:i+k_length])
        return k_mers
    
#Need to generate random kmers that aren't in the chromosome so 
#that I can actually test for false positives in the bloom filter later. 
def GenerateRandomKmers(k_mers, num_of_k_mers, length):
    not_in_kmers = []
    nucls = ["A", "C", "T", "G"]
    while len(not_in_kmers) < num_of_k_mers:
        new_kmer = ''
        for k in list(np.random.choice(nucls, length)):
            new_kmer += k
        if new_kmer in k_mers:
            break
        else:
            not_in_kmers.append(new_kmer)
    return not_in_kmers


##Just a quick check to see if things are adding up.
def Check_If_Bloom_Works(k_mers, not_in_kmers, bloomfilter):
    np.random.shuffle(k_mers)
    np.random.shuffle(not_in_kmers)
    test_mers = k_mers[:1000] + not_in_kmers[:75]
    np.random.shuffle(test_mers)
    false_pos = 0
    true_pos = 0
    true_false = 0
    for i in test_mers:
        ret_var = bloomfilter.check(i)
        if ret_var == True and i in not_in_kmers:
            false_pos += 1
        elif ret_var == False and i in not_in_kmers:
            true_false += 1
        elif ret_var == True and i in k_mers:
            true_pos += 1
    print("True positives were: {}".format(true_pos))
    print("True Falses were: {}".format(true_false))
    print("False_positives were: {}".format(false_pos))
            


    
    

    
    
