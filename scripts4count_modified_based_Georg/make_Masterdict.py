__author__ = 'georg.michlits'

import sys
import os

print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)
print(sys.argv[1]) # input
print(sys.argv[2]) # guide file
print(sys.argv[3]) # experiment index
print(sys.argv[4]) # ouptput

#bowtie_result_file = open(sys.argv[1], 'r')
input_out0_file = open(sys.argv[1], 'r')

guidefile_name = sys.argv[2]
indexfilename = sys.argv[3]
sgALLfile = open(guidefile_name,'r')
Indexfile = open(indexfilename,'r')

output_filename = sys.argv[4]
output_2_MD = open(output_filename,'w')
output_1_kicked_f = open(output_filename + '_trashed_.txt','w')

# mismatches (mm) for sample barcodes 
index_mm_allowed = 0
# missmatches for guides
guide_mm_allowed = 1
#index_mm_allowed = int(input('enter number of index MM allowed (0): '))
#guide_mm_allowed = int(input('enter number of guide MM allowed (1): '))




#create index dictonaries that contain information on indices

def MMmax(dna1,dna2,maxMM):
    c=0
    for i in range(0,len(dna1)):
        if not dna1[i] == dna2[i]:
            c = c + 1
            if c > maxMM:
                return False
                break
    return True

Indexdict = {}
for line in Indexfile:
    element = line.rstrip('\n').split(',')
    Exp_name = element[1]
    Index_n = element[5]
    Indexdict[Index_n] = Exp_name
print('Index dictionaries done')
Indexfile.close()

#from file containing all guides create a dict that later counts in experiments - clones - count
gDNA = {}
readout = {}
for line in sgALLfile:
    element = line.rstrip('\n').split('\t')
    gDNAname = (element[0])
    guide20nt = element[1]
    gDNA[gDNAname] = guide20nt
    readout[gDNAname] = {}
sgALLfile.close()

good_index_good_guide = 0
count_kicked_lines = 0
bad_index_kicked = 0
guide_mm_kicked = 0
total_lines = 0

ruled_out_index = 0
print('reading NGS.sam file')
print('generating MD_dict')
for line in input_out0_file:
    if len(line.rstrip('\n').split('\t')) == 6:
        total_lines = total_lines + 1
        if total_lines%10000000 == 0:
            print(str(total_lines/1000000) + ' million')
        #print(total_lines)
        element = line.rstrip('\n').split('\t')
        gDNAname = element[4]
        guide_mm = int(element[5])
        bc2_index = element[2][5:11]
        barcode = element[3][5:15]
        test_if_index_ok = 'n'
        test_if_guide_mm_ok = 'n'
        if guide_mm <= guide_mm_allowed:
            test_if_guide_mm_ok = 'y'
        for index in Indexdict:
            if test_if_guide_mm_ok == 'y' and MMmax(bc2_index,index,index_mm_allowed):
                test_if_index_ok = 'y'
                good_index_good_guide = good_index_good_guide + 1
                exp_name = Indexdict[index]
                #guide20nt = gDNA[gDNAname]
                #count it into the readoutdictionary
                #readout = {'gDNAname:{index:{barcode:count}}}
                if not exp_name in readout[gDNAname]:
                    readout[gDNAname][exp_name] = {}
                if not barcode in readout[gDNAname][exp_name]:
                    readout[gDNAname][exp_name][barcode] = 0
                readout[gDNAname][exp_name][barcode] = readout[gDNAname][exp_name][barcode] + 1
                break
        if test_if_index_ok == 'n' or test_if_guide_mm_ok == 'n':
            output_1_kicked_f.write(line)
            count_kicked_lines = count_kicked_lines + 1
        if test_if_index_ok == 'n':
            bad_index_kicked += 1
        if test_if_guide_mm_ok == 'n':
            guide_mm_kicked += 1
input_out0_file.close()

print('total lines: ' + str(total_lines))
print('good lines (index and guide match <= max MM): ' + str(good_index_good_guide))
print('kicked lines: ' + str(count_kicked_lines))
print('kicked (bad index): ' + str(bad_index_kicked))
print('kicked (guide MM): ' + str(guide_mm_kicked))

# write output file (the saved file can later be used to read in the dict again)

print('writing Masterdict')
for gDNAname in sorted(readout):
    for exp_name in sorted(readout[gDNAname]):
        for i in sorted(readout[gDNAname][exp_name].items(), key=lambda t: t[1], reverse=True):
            barcode = i[0]
            c = i[1]
            output_2_MD.write(gDNAname + '\t' + exp_name + '\t' + barcode + '\t' + str(c) + '\n')

print('writing Masterdict finished')
output_2_MD.close()