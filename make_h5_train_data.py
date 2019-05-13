
# coding: utf-8

# Reading VCF file and preparing data for VCF Challenge

# In[97]:


# load modlues
import numpy as np

import re
import csv
import h5py


# In[95]:


def splitRa(ra):
    # Helper function to split RA field into forward_and_reverse_read_count, mapping_qualities, pile_up grid
    # Input RA filed as list as read by vcf parser
    # Outputs:
    #   forward_reverse_reads: 3 x 2 numpy array
    #   mapping_qs: 3 x 6 numpy array
    #   pile_up: 3 x 61 x 13 numpy array
    # If only 2 bins present will copy bin1 to bin and much 2 to 3 to maintain shape
    
    # collapse to string remove "|"
    str_ra = ra
    str_ra = re.sub(r",", "*", str_ra)
    str_ra = re.sub(r"\|", "*", str_ra)
    str_ra = re.sub(r"\:", "*", str_ra)
    str_ra = re.sub(r"\{", "", str_ra)

    split_ra = str_ra.split("}")  # now split bins on "{"
    split_ra = split_ra[0:(len(split_ra)-1)]

    # split back into lists
    ra = [x.split("*") for x in split_ra]
    
    num_bins = len(ra)  # get number of bins

    forward_reverse_reads = np.zeros([3, 2])
    mapping_qs = np.zeros([3, 6])
    pile_up = np.zeros([3, 61, 13])

    for i in range(len(ra)):
        rf = ra[i]
        rf = [float(x) for x in rf] # covnert to float
        forward_reverse_reads[i, 0:2] = rf[0:2]
        mapping_qs[i, 0:7] = rf[2:8]
        for j in range(13):
            pile_up[i, 0:62, j] = rf[(8+(j*61)):(69+(j*61))]
            
    # if only 2 bins copy first bin into row2 and push 2 to 3
    if num_bins == 2:
        forward_reverse_reads[2,] = forward_reverse_reads[1,]
        forward_reverse_reads[1,] = forward_reverse_reads[0,]
        mapping_qs[2,] = mapping_qs[1,]
        mapping_qs[1,] = mapping_qs[0,]
        pile_up[2,] = pile_up[1,]
        pile_up[1,] = pile_up[0,]

    return forward_reverse_reads, mapping_qs, pile_up

# FUNCTION to READ FILED FROM ENTIRE DATA
# Read Fields into seprate numpy arrays
def readVcfFields(file):
    crf = file[0]
    frf = file[1]
    gc = file[2]
    frr, mqs, pu = splitRa(''.join(file[3]))

    return crf, frf, gc ,frr ,mqs, pu


# In[179]:


# TEST FALSE POSITIVE ===============================================
#set dimensions / create hdf5 file
nrows_fp = 41255

h5f = h5py.File("./test_data_LC_fp.h5", 'w')

set_crf = h5f.create_dataset('CRF', (nrows_fp,), dtype='float')
set_frf = h5f.create_dataset('FRF', (nrows_fp,), dtype='float')
set_gc = h5f.create_dataset('GC', (nrows_fp,), dtype='float')
set_frr = h5f.create_dataset('FRR', (nrows_fp, 3, 2), dtype='float')
set_mq = h5f.create_dataset('MQ', (nrows_fp, 3, 6), dtype='float')
set_pu = h5f.create_dataset('PU', (nrows_fp, 3, 61, 13), dtype='i')
h5f.create_dataset('LABELS', data=np.ones(nrows_fp), dtype='i')

in_file = "./data/test/octopus.hackathon.1000G.LC.CSR.fp.test.tab"

# fill with pileup tensor
i = 0
with open(in_file, 'r') as file:
    for line in csv.reader(file, delimiter="\t"):
        crf, frf, gc ,frr ,mq, pu = readVcfFields(line)
        set_crf[i,] = float(crf)
        set_frf[i,] = float(frf)
        set_gc[i,] = float(gc)
        set_frr[i,] = frr
        set_mq[i,] = mq
        set_pu[i,] = pu
        i += 1
        if i >= nrows_fp:
            break

h5f.close()


# In[180]:


# TEST TRUE POSITIVE ===============================================
#set dimensions / create hdf5 file
nrows_fp = 318865

h5f = h5py.File("./test_data_LC_tp.h5", 'w')

set_crf = h5f.create_dataset('CRF', (nrows_fp,), dtype='float')
set_frf = h5f.create_dataset('FRF', (nrows_fp,), dtype='float')
set_gc = h5f.create_dataset('GC', (nrows_fp,), dtype='float')
set_frr = h5f.create_dataset('FRR', (nrows_fp, 3, 2), dtype='float')
set_mq = h5f.create_dataset('MQ', (nrows_fp, 3, 6), dtype='float')
set_pu = h5f.create_dataset('PU', (nrows_fp, 3, 61, 13), dtype='i')
h5f.create_dataset('LABELS', data=np.ones(nrows_fp), dtype='i')

in_file = "./data/test/octopus.hackathon.1000G.LC.CSR.tp.test.tab"

# fill with pileup tensor
i = 0
with open(in_file, 'r') as file:
    for line in csv.reader(file, delimiter="\t"):
        crf, frf, gc ,frr ,mq, pu = readVcfFields(line)
        set_crf[i,] = float(crf)
        set_frf[i,] = float(frf)
        set_gc[i,] = float(gc)
        set_frr[i,] = frr
        set_mq[i,] = mq
        set_pu[i,] = pu
        i += 1
        if i >= nrows_fp:
            break

h5f.close()

