#!/usr/bin/env python
# coding: utf-8

# # Complete sequencing of barcoded pol mutants
# 
# Allie Stanton
# 
# Data is from a 500 cycle MiSeq run with 15 million reads started on April 5, 2019. The input material was from 293Ts transfected with 0.1 MOI of lentivirus with a 15N barcode (plasmid 2626) with different single AA polymerase mutants. The cells were flow-sorted on BFP to isolate 10k infected cells from each sample. The 10k cells from each sample were split into 10 pools of 1000 cells after gDNA extraction, such that each sample has 10 indices associated with it.

# ## Reading in files
# 
# First I'll create a dictionary with the sample names (mutants) as keys and read 1 and read 2 file names as values. Then I'll define the functions that read in the data. I'll use the WT sample as an example to go through each of the functions once and demonstrate their purpose before processing all of the samples. Note that each sample has 10 indices associated with it so I'll read them in on a loop, keeping the indices separate for now.

# In[1]:


from Bio import SeqIO
import itertools
import os
import random
import matplotlib.pyplot as plt
import numpy as np
import scipy

files = {}
seqdir = os.path.join(os.getcwd(), "sequences")
for directory in os.listdir(seqdir):
    if not directory.startswith("."):
        longdir = os.path.join(seqdir, directory)
        for file in os.listdir(longdir):
            filename = os.fsdecode(file)
            if file.endswith("R1_001.fastq"):
                read1_file = os.path.join(longdir, filename)
            elif file.endswith("R2_001.fastq"):
                read2_file = os.path.join(longdir, filename)
        files[directory] = (read1_file, read2_file)


# In[2]:


def read1_in(file):
######################################################################################
#                                                                                    #
# Reads in the read 1 sequence data and stores sequence and barcode data separately. #    
#                                                                                    #
# Inputs:                                                                            #
#      file: full path of file where read 1 data is stored in fastq format.          #
# Outputs:                                                                           #
#      barcodes: dictionary with record ids as keys and a string containing the 15N  #
#        barcode as values.                                                          #          
#      read1_seqs: dictionary with record ids as keys and a string containing the    #
#        first 100 bases of the corresponding read 1 after the barcode as values.    #
#                                                                                    #
######################################################################################

    barcodes = {}
    read1_seqs = {}
    for record in SeqIO.parse(file, "fastq"):
        # Make sure the adapter is present and sequence is correct length
        if str(record.seq[9:24]) == 'GAGCTGTACAAATAA' and str(record.seq[39:44]) == 'GGATC' and len(str(record.seq)) >= 248:            
            q_scores = record.letter_annotations["phred_quality"]
            q_scores_read = q_scores[39:248]
            read = list(str(record.seq[39:248]))
            
            bc = record.seq[24:39]

            read_qfilter = ''
            n_count = 0
            for i in range(len(read)):
                if q_scores_read[i] < 30:
                    read_qfilter += 'N'
                    n_count += 1
                else:
                    read_qfilter += read[i]
            
            if n_count <= 20:
                read1_seqs[record.id] = read_qfilter
                barcodes[record.id] = str(bc)
    
    return barcodes, read1_seqs


# In[3]:


WT_barcodes = []
WT_read1_seqs = []
nreads = 0
for i in range(1, 11):
    folder = 'WT' + str(i)
    file_subset = files[folder]
    barcodes, read1_seqs = read1_in(file_subset[0])
    
    nreads += len(read1_seqs)
    
    WT_barcodes.append(barcodes)
    WT_read1_seqs.append(read1_seqs)
    
print(nreads, 'valid read 1 reads')


# In[4]:


def read2_in(file):
######################################################################################
#                                                                                    #
# Reads in the read 2 sequence data and stores sequence and barcode data separately. #    
#                                                                                    #
# Inputs:                                                                            #
#      file: full path of file where read 2 data is stored in fastq format.          #
# Outputs:                                                                           #
#      read2_seqs: dictionary with record ids as keys and a string containing the    #
#        first 100 bases of the corresponding read 2 as values.                      #
#                                                                                    #
######################################################################################

    read2_seqs = {}
    for record in SeqIO.parse(file, "fastq"):
        # Make sure the adapter is present and sequence is the correct length
        if str(record.seq[3:21]) == 'ACCCAAAGTAGTCGGTTC' and len(str(record.seq)) >= 248:
            q_scores = record.letter_annotations["phred_quality"]
            q_scores_read = q_scores[21:248]
            
            read = list(str(record.seq[21:248]))
            
            read_qfilter = ''
            n_count = 0
            for i in range(len(read)):
                if q_scores_read[i] < 30:
                    read_qfilter += 'N'
                    n_count += 1
                else:
                    read_qfilter += read[i]
            
            if n_count <= 20:
                read2_seqs[record.id] = read_qfilter
    
    return read2_seqs


# In[5]:


WT_read2_seqs = []
nreads = 0
for i in range(1, 11):
    folder = 'WT' + str(i)
    file_subset = files[folder]
    read2_seqs = read2_in(file_subset[1])
    
    nreads += len(read2_seqs)

    WT_read2_seqs.append(read2_seqs)
    
print(nreads, 'valid read 2 reads')


# ## Grouping read pairs into unique infection events by barcode
# Each half of the read pairs is now associated with a unique record id that corresponds to a single molecule. Moreover, each barcode is now associated with one or more record ids. Most barcodes probably have more than one read pair associated with them as the DNA molecules from each infection are amplified by PCR. We now want to take each unique barcode and, using the record id as an intermediate, find all of the read pairs that are associated with that barcode and group them together. Meanwhile, since we filtered out some reads in the previous step that had adaptor or length issues, we want to make sure that we only keep read pairs where both halves of the pair exist. 
# 
# First we want to focus on just all unique barcodes and make groups of reads associated with identical barcodes, which should correspond to one integration associated with one infection event (more on this assumption later). Then we will collapse these groups if the barcodes are very similar. Finally, we will (somewhat heuristically) remove any abnormally large or small groups.

# In[6]:


def find_barcode_groups(barcodes, read1_seqs, read2_seqs):
######################################################################################
#                                                                                    #
# Groups together read pairs that have the same barcode.                             #    
#                                                                                    #
# Inputs:                                                                            #
#      barcodes: dictionary with record ids as keys and a string containing the 15N  #
#        barcode as values.                                                          #  
#      read1_seqs: dictionary with record ids as keys and a string containing the    #
#        first 100 bases of the corresponding read 1 after the barcode as values.    #
#      read2_seqs: dictionary with record ids as keys and a string containing the    #
#        first 100 bases of the corresponding read 2 as values.                      #
# Outputs:                                                                           #
#      barcode_groups: dictionary with 15N barcode (type: str) as keys and           #
#        corresponding read pairs (type: tuple of strings) as values. Each key may   #
#        have multiple values.                                                       #
#      unique_valid: list of unique barcodes (list of strings).                      #
#                                                                                    #
######################################################################################

    # intersect the record ids to find record ids that have both read 1 and read 2
    shared_ids = list(read1_seqs.keys() & read2_seqs.keys())

    barcode_groups = {}
    for idx in shared_ids:
        readpair = (read1_seqs[idx], read2_seqs[idx])
        bc = barcodes[idx]
        
        # create a new sublist for each unique barcode and fill it with read pair tuple
        barcode_groups[bc] = barcode_groups.get(bc, [])
        barcode_groups[bc].append(readpair)
    
    # get a list of unique barcodes
    unique_valid = list(barcode_groups.keys())
    
    return barcode_groups, unique_valid


# In[7]:


WT_barcode_groups = []
WT_unique_valid = []
ngroups = 0
for i in range(0, 10):
    barcode_groups, unique_valid = find_barcode_groups(WT_barcodes[i], WT_read1_seqs[i], WT_read2_seqs[i])
    ngroups += len(barcode_groups)
    WT_barcode_groups.append(barcode_groups)
    WT_unique_valid.append(unique_valid)

print(ngroups, "groups of valid read pairs where each group has a unique barcode")


# ### Hamming distance
# Hamming distance is a measure of the difference between two strings of equal length. It is equal to the number of positions at which the corresponding symbols are different. In this implementation, we calculate the Hamming distance by counting the number of unequal positions between the two strings. An example Hamming distance is given between two barcodes drawn from the WT sample.

# In[8]:


def hamming(s1, s2):
######################################################################################
#                                                                                    #
# Takes two sequences and computes the Hamming distance between them.                #  
#                                                                                    #
# Inputs:                                                                            #
#      s1, s2: two DNA sequences of equal length                                     #
# Outputs:                                                                           #
#      dist: an integer representing the hamming distance                            #
#                                                                                    #
######################################################################################
    
    if isinstance(s1, str) == 0:
        s1 = str(s1)
    if isinstance(s2, str) == 0:
        s2 = str(s2)
    
    s1 = list(s1)
    s2 = list(s2)

    dist = len([i for i, j in zip(s1, s2) if i != j and i != 'N' and j != 'N'])
    
    return dist


# In[9]:


print(WT_unique_valid[0][0])
print(WT_unique_valid[0][1])
print("Hamming distance:", hamming(WT_unique_valid[0][0], WT_unique_valid[0][1]))


# ### Probabilistic interlude #1
# We are about to use the Hamming distance metric to collapse barcode groups where the barcodes are very similar—specifically, when the distance between two barcodes is one. The motivation for this is that it seems relatively rare that two random 15mers would be so close as to differ by only one base, and PCR errors may be responsible for modifying barcodes that in fact originated from the same infection event.
# 
# But how likely exactly is it that without errors, no pair of (presumably) random barcodes in a group of several thousand barcodes are exactly one base apart? Well, we can easily find this mathematically: given an equal probability of mutation from or to any base, the probability that two given random 15mers differ by exactly one base is 
# 
# $15\times 0.25^{14}\times 0.75 = 4.191\times 10^{-8}$
# 
# Extend this to 1000 barcodes (the theoretical number of barcodes in each index) and the probability of having all barcodes be at least 2 bases apart from all other barcodes is
# 
# $(1 - (15\times (0.25^{14}\times 0.75))^{1000} = 0.999958$
# 
# It therefore seems safe to assume any barcodes that are 1 base apart are a result of PCR errors and did not truly originate from different events.
# 
# Note that in the above hamming_collapse function, I collapsed the barcodes into groups of a particular geometry. By this I mean that if a group is depicted as a graph where barcodes are nodes and edges exist where two nodes are at a Hamming distance of one, the graph has the following properties: 
# 
# 1. Each node shares an edge with at least one other node
# 
# 2. Every pair of nodes has a path between them (the graph is connected)
# 
# 3. Each pair of nodes does not necessarily have an edge between them (the graph is not necessarily complete)
# 
# In other words, if a group contains barcodes A, B, and C, A and B may be different by one base and B and C may be different by one base, but A and C can be different by either one or two bases. This can theoretically result in an arbitrarily long chain of barcodes. I chose to ignore this (instead removing large groups at the end) as I couldnt come up with a better solution. 
# 

# In[10]:


def hamming_collapse(unique_valid): 
######################################################################################
#                                                                                    #
# Collapses barcodes into groups where barcodes within the group differ by one base  #
#   (Hamming distance of 1) and removes groups that contain more than 3 barcodes to  #
#   remove abnormally complex groups.                                                #
#                                                                                    #
# Inputs:                                                                            #
#      unique_valid: list of unique barcodes (list of strings).                      #
# Outputs:                                                                           #
#      hamm1_reduced: list of sublists where each sublist contains all of the tuples #
#        of read 1 and read 2 pairs that belong in one barcode group. If the barcode #
#        group has more than one barcode, the sublist will contain as many lists as  #
#        barcodes in the group.                                                      #
#                                                                                    #
######################################################################################

    unique_hamm1 = dict((k, 0) for k in unique_valid)
    unique_hamm1[unique_valid[0]] = 1

    counter = 0
    for i in unique_valid:
        counter += 1
        if unique_hamm1[i] == 0:
            unique_hamm1[i] = max(unique_hamm1.values()) + 1
        for j in unique_valid[counter: ]:
            if hamming(i, j) <= 1:
                unique_hamm1[j] = unique_hamm1[i]

    hamm1_groups = {}
    for k, v in unique_hamm1.items():
        hamm1_groups[v] = hamm1_groups.get(v, [])
        hamm1_groups[v].append(k)
    
    hamm1_reduced = []
    for v in hamm1_groups.values():
        if len(v) <= 6:
            hamm1_reduced.append(v)
    
    return hamm1_reduced


# In[11]:


WT_hamm1_reduced = []
ngroups_collapsed = 0
for i in range(0, 10):
    hamm1_reduced = hamming_collapse(WT_unique_valid[i])
    ngroups_collapsed += len(hamm1_reduced)
    WT_hamm1_reduced.append(hamm1_reduced)

print(ngroups_collapsed, "barcodes after collapsing within a Hamming distance of 1 and removing groups of large size")


# ## Flattening and quality control
# Now that we have the barcode groups, all that remains are a few quality control and processing steps. First we need to flatten the sublists in each barcode group (to remove the boundary between two barcodes within the same group) and remove barcode groups with a very small number of read pairs. Then we will concatenate the read 1 and read 2 sequences into one string, which will make them slightly easier to work with later. Then we will make sure that each sequence within a group is identical (more on this in a bit).
# 
# First, let's take a look at the number of read pairs per barcode.

# In[25]:


reads_per_barcode = []
for i in range(0,10):
    for group in hamm1_reduced:
        seqs_hamm = []
        for i in range(len(group)):
            seqs_hamm.append(barcode_groups[group[i]])
        seqs_hamm_flat = [seq for sublist in seqs_hamm for seq in sublist]
        reads_per_barcode.append(len(seqs_hamm_flat))

_, _, patch = plt.hist(reads_per_barcode, 100, density=1, histtype = 'step', facecolor='blue')

plt.xlabel('Number of reads')
plt.ylabel('Probability')
plt.title('Histogram of Number of Reads per Barcode')
plt.show()


# As we can see, there are a lot of barcodes with very few (~1) reads. This is probably some sort of artifact of PCR or sequencing error. Let's see what happens if we plot only the barcodes with at least 10 reads:

# In[28]:


reads_per_barcode = []
for i in range(0,10):
    for group in hamm1_reduced:
        seqs_hamm = []
        for i in range(len(group)):
            seqs_hamm.append(barcode_groups[group[i]])
        seqs_hamm_flat = [seq for sublist in seqs_hamm for seq in sublist]
        
        if len(seqs_hamm_flat) >= 10:
            reads_per_barcode.append(len(seqs_hamm_flat))

_, _, patch = plt.hist(reads_per_barcode, 100, density=1, histtype = 'step', facecolor='blue')

plt.xlabel('Number of reads')
plt.ylabel('Probability')
plt.title('Histogram of Number of Reads per Barcode')
plt.show()


# This looks a lot more reasonable, and it's also no longer zero-centered (or 10-centered, as you will), which indicates that the "real" barcodes with >10 reads are saturated. In other words, if we did more sequencing on this data, we would probably not recover any new real barcodes. Now let's use the cutoff of 10 reads per barcode.

# In[30]:


def combine_seqs(hamm1_reduced, barcode_groups):
######################################################################################
#                                                                                    #
# A housekeeping function that flattens the barcode groups into a list of sublists   # 
#   where each sublist contains the reads for one barcode group, removes barcode     #
#   groups that don't have at least two read pairs, and finally concatenates the     #
#   read pair into one string that contains both read 1 and read 2.                  #
#                                                                                    #
# Inputs:                                                                            #
#      hamm1_reduced: list of sublists where each sublist contains all of the tuples #
#        of read 1 and read 2 pairs that belong in one barcode group. If the barcode #
#        group has more than one barcode, the sublist will contain as many lists as  #
#        barcodes in the group.                                                      #
#      barcode_groups: dictionary with 15N barcode (type: str) as keys and           #
#        corresponding read pairs (type: tuple of strings) as values. Each key may   #
#        have multiple values.                                                       #
# Outputs:                                                                           #
#      seqs_concat: a list of sublists where each sublist contains all of the        #
#        sequences in a barcode group (list of flattened sublists of strings). Each  #
#        string corresponds to one read pair.                                        #
#                                                                                    #
######################################################################################

    seqs_groups = []
    for group in hamm1_reduced:
        seqs_hamm = []
        for i in range(len(group)):
            seqs_hamm.append(barcode_groups[group[i]])
        seqs_hamm_flat = [seq for sublist in seqs_hamm for seq in sublist]
        if len(seqs_hamm_flat) >= 10:
            seqs_groups.append(seqs_hamm_flat)
    
    seqs_concat = []
    for group in seqs_groups:
        seqs_concat_grp = []
        for pair in group:
            seqs_concat_grp.append(pair[0] + pair[1])
        seqs_concat.append(seqs_concat_grp)
        
    return seqs_concat


# In[32]:


WT_seqs_concat_nested = []
for i in range(0, 10):
    seqs_concat = combine_seqs(WT_hamm1_reduced[i], WT_barcode_groups[i])
    WT_seqs_concat_nested.append(seqs_concat)

WT_seqs_concat = [seq for sublist in WT_seqs_concat_nested for seq in sublist]

print(len(WT_seqs_concat), "barcode groups with ten or more read pairs")


# ### Managing internal disagreement within barcode groups
# So far we have created groups of reads that either came from the same barcode or came from barcodes that were only one base apart. We assume that each of these groups of reads corresponds to a single unique virus infecting a single cell, which should mean that all of the sequences within a barcode group are identical save for errors introduced in PCR. Note that this assumption only holds if there were the correct number of viruses added to the cells such that 
# 
# 1. Each cell is only infected by one virus, and
# 
# 2. No two viruses have the same barcode.
# 
# 1 is mostly taken care of by the low MOI used to infect (MOI of ~0.1 ensures <1% multiple infection). 2 is trickier and it relies on the assumption that the library of viruses is sufficiently diverse to favor unique barcodes. My calculations indicated that around 57,000 unique transformants went into the maxiprep that produced the plasmid containing the barcode. Assuming that this complexity was preserved through maxiprep and transfection, at an MOI of 0.1 each well containing 2,000,000 cells was infected with 100,000 viruses. The number of times each barcode is drawn with 100,000 trials can be modeled with a multinomial distribution where the probability of each barcode being drawn is 1/57000. Let's model this in the code block below.

# In[33]:


mv = np.random.multinomial(100000, [1/57000]*57000, 100)
total_viruses = np.sum(mv >= 1, 1)
total_multiples = np.sum(mv >= 2, 1)
percent_multiples = np.average(np.divide(total_multiples, total_viruses))
print(round(100*percent_multiples, 2), "% of all barcodes correspond to non-unique barcodes")


# However, we aren't taking all of the cells in the original sample, only 10,000 of them. This introduces some narrowing of barcodes and already starts to limit the number of barcode collisions. We can model this with a multivariate hypergeometric distribution, which is the generalization of the hypergeometric distribution. It goes like this: say that there are $K_{i}$ balls of color $i$ in the urn, and you draw $n$ balls without replacement. Then the number of balls of each color in the sample $(k_{1}, k_{2}, ..., k_{c})$ has the multivariate hypergeometric distribution. In python we define this distribution in terms of the hypergeometric distribution.

# In[34]:


def multi_hypergeom(m, colors):
    remaining = np.cumsum(colors[::-1])[::-1]
    result = np.zeros(len(colors), dtype=np.int)
    for i in range(len(colors)-1):
        if m < 1:
            break
        result[i] = np.random.hypergeometric(colors[i], remaining[i+1], m)
        m -= result[i]
    result[-1] = m
    return result


# Now the colors are barcodes and the balls are viruses, and the number of viruses with each barcode is taken from the multinomial sampling done above. First, we do 10,000 draws to model the initial flow sorting of 10,000 cells. Then we do 1000 draws on the result of the original 10,000 draws to model partitioning the cells out into ten different indices.

# In[38]:


virus_counts_unstripped = list(mv[0])
virus_counts = []
for i in virus_counts_unstripped:
    if i != 0:
        virus_counts.append(i)

# Sample 10,000 cells from flow sorting
virus_counts = np.array(virus_counts)
draws_sort = 10000
counts_sorted = multi_hypergeom(draws_sort, virus_counts)

# Sample 1000 from the 10,000 to model partitioning into indices
draws_partition = 1000
counts_partitioned = multi_hypergeom(draws_partition, counts_sorted)

total_viruses = np.sum(counts_partitioned)
total_multiples = np.sum(counts_partitioned >= 2)
percent_multiples = total_multiples / total_viruses
print(round(100*percent_multiples, 2), "% of all barcodes correspond to non-unique barcodes")


# We now know that very few barcodes will be examples of two identical viruses infecting different cells. However, some of the reads contained within a group of reads corresponding to a single barcode still have internal disagreements. This is most likely due to PCR errors. We therefore need to correct these errors where we can by identifying the consensus sequence across reads belonging to one barcode, discarding barcodes (and corresponding reads) when no consensus can be found.
# 
# In order to do that, I am going to take the average base at each position and make that the consensus base at said position. If any position does not have at least 50% of the reads as one base, then that entire group will be discarded. In order to find the average base we must first convert sequences to one-hot encoding, and then find the consensus sequence by base averaging.

# In[39]:


def sequences_to_one_hot(sequences, chars = 'ACGTN'):
    
    seqlen = len(sequences[0])
    char_to_int = dict((c, i) for i, c in enumerate(chars))

    one_hot_encoded = []
    for seq in sequences:
        onehot_seq = []
        
        integer_encoded = [char_to_int[base] for base in seq]
        for value in integer_encoded:
            letter = [0 for _ in range(len(chars))]
            letter[value] = 1
            onehot_seq.append(letter)
        one_hot_encoded.append(onehot_seq)
        
    one_hot_encoded = np.array(one_hot_encoded)
    
    return one_hot_encoded


# In[40]:


def find_multiple_alignment(seqs_concat):
    
    true_seqs = []
    chars = 'ACGTN'
    for group in seqs_concat:
        onehot_group = sequences_to_one_hot(group, chars)
        onehot_group = np.delete(onehot_group, 4, 2)
        
        base_averages = np.average(onehot_group, axis = 0)
        top_base = np.argmax(base_averages, axis = 1)
        top_base_acc = np.amax(base_averages, axis = 1)

        int_to_char = dict((i, c) for i, c in enumerate(chars))
        
        if np.amin(top_base_acc) >= 0.50:
            consensus_bases = np.ndarray.tolist(top_base)
            consensus_bases = [int_to_char[base] for base in consensus_bases]
            true_seqs.append(''.join(consensus_bases))
    
    return true_seqs


# In[28]:


WT_true_seqs = find_multiple_alignment(WT_seqs_concat)
print(len(WT_seqs_concat) - len(WT_true_seqs), "barcode groups with mismatches between reads")
print(len(WT_true_seqs), "barcode groups without mismatches between reads")


# ## Error rate calculation
# Now we have sequences that theoretically correspond to single infection events. All that's left to calculate the error rate by comparing the sequences to the CMV IRES (the region that was sequenced).

# In[29]:


IRES_read1_fasta = SeqIO.read("IRES_read1.fasta", "fasta")
IRES_read1 = IRES_read1_fasta.seq[0:209]

IRES_read2_fasta = SeqIO.read("IRES_read2.fasta", "fasta")
IRES_read2 = IRES_read2_fasta.seq[21:248]

IRES = str(IRES_read1) + str(IRES_read2)


# In[30]:


def calc_error_rate(true_seqs, reference):
######################################################################################
#                                                                                    #
# Finds errors between sequences and a reference and calculates the total error      #
#   rate.                                                                            #
#                                                                                    #
# Inputs:                                                                            #
#      true_seqs: a list of sequences where each sequence corresponds one-to-one     #
#        with a barcode group. List of strings.                                      #
#      reference: a string of the same length of each element of true_seqs           #
#        representing the reference sequence.                                        #
# Outputs:                                                                           #
#      bases_per_error: float representing the error rate in bases per error.        #
#                                                                                    #
######################################################################################

    nerror_list = []
    seqs_tested = 0
    for seq in true_seqs:
        errors = hamming(seq, reference)
        if errors <= 2:
            nerror_list.append(errors)
            seqs_tested += 1
        
    nerrors = sum(nerror_list)
    bases_tested = len(reference) * seqs_tested
    if sum(nerror_list) == 0:
        bases_per_error = float('nan')
    else:
        bases_per_error =  nerrors / bases_tested
    
    return bases_per_error, nerrors, bases_tested, seqs_tested


# In[31]:


WT_error_rate, _, _, _ = calc_error_rate(WT_true_seqs, IRES)
print("Error rate: 1 error per", int(1/WT_error_rate), "bases")


# In[32]:


barcode_numbers = {}
error_rates = {}
samples = ['WT', 'A272S', 'K65R', 'K70E', 'L74I', 'L74V', 'M184I', 'Q151M']

for sample in samples:
    print("Calculating for sample", sample)
    
    s_seqs_concat_nested = []
    for i in range(1, 11):
        folder = sample + str(i)
        file_subset = files[folder]
        
        barcodes, read1_seqs = read1_in(file_subset[0])
        read2_seqs = read2_in(file_subset[1])
        barcode_groups, unique_valid = find_barcode_groups(barcodes, read1_seqs, read2_seqs)
        
        hamm1_reduced = hamming_collapse(unique_valid)
        
        seqs_concat = combine_seqs(hamm1_reduced, barcode_groups)
        s_seqs_concat_nested.append(seqs_concat)
    
    s_seqs_concat = [seq for sublist in s_seqs_concat_nested for seq in sublist]  
    s_true_seqs = find_multiple_alignment(s_seqs_concat)
    
    s_error_rate, s_nerrors, s_bases_tested, s_seqs_tested = calc_error_rate(s_true_seqs, IRES)
    
    barcode_numbers[sample] = s_seqs_tested
    error_rates[sample] = (s_error_rate, s_nerrors, s_bases_tested)
        
print("Done")


# In[33]:


print(barcode_numbers)


# In[34]:


print(error_rates)

