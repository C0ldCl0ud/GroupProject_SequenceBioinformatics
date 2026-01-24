#!/usr/bin/env python3

#==============================================================================
# EVALUATION SCRIPT - Group Project SeqBI WS25
#
# execute via shell script: ./evaluate.sh
#
#
# IMPORTANT:
# copy this script to each ".../results/"
# and run from there with each dataset/tool in
# in DATASETS/TOOLS uncommented, which is
# present in that specific results folder.
#
# ADD ALL KNOWN TOOLS WHICH PRODUCE .fasta
# INSTEAD OF .fa FILES TO THE LIST BELOW
#==============================================================================

import csv


#=======================================#
# MAKE CHANGES TO THIS 3 variable LISTS #
#=======================================#

DATASETS = [
    'cheese',
#    'marine',
#    'activated_sludge',
#    'human_gut_1',
#    'human_gut_2'
    ]

TOOLS = [
    'metabat2',
    'maxbin2',
    'concoct',
#    'vamb',
#    'metadecoder',
#    'binny',
#    'metabinner',
    'semibin2',
#    'comebin'
    ]

# add all known tools, that name their binnings with the ending .fasta instead of .fa
NOT_FA_BUT_FASTA = [
    'maxbin2',
    ]


#==================================================#
# THESE global parameters DON'T HAVE TO BE CHANGED #
#==================================================#

SAMPLING_TYPE = ['single', 'multi', 'coassembly']
ASSEMBLY_TYPE = ['short', 'long', 'hybrid']
# S00 .. S07
SAMPLES = [f'S0{i}' for i in range(8)]
# sampling assembly combinations (like in the heatmap plot)
MODES = [
    'short_co',
    'short_single',
    'short_multi',
    'long_single',
    'long_multi',
    'hybrid_single',
    'hybrid_multi'
    ]
# bin quality scheme
BIN_QUALITIES = ['HQ','NC','MQ', 'total_bin_count']



#=============#
# PREPARATION #
#=============#


# reads in the list of available paths to evaluate (created by evaluate.sh)
with open('paths.eval', newline='', encoding='utf-8') as paths_f:
    RESULT_PATHS = paths_f.read().splitlines()

# creates a nested (4-dim) dictionary hierarchically:
    # 1. dataset
    # 2. tool
    # 3. sampling-type assembly-type combination
    # 4. bin quality
eval_bin_counts_dict = dict()

for d in DATASETS:
    eval_bin_counts_dict[d] = dict()
    for t in TOOLS:
        eval_bin_counts_dict[d][t] = dict()
        for m in MODES:
            eval_bin_counts_dict[d][t][m] = dict()
            for bq in BIN_QUALITIES:
                eval_bin_counts_dict[d][t][m][bq] = 0




#======================================#
# EVAL QUALITY OF BINNINGS - FIND MAGs #
#======================================#

MAGs = list()


print('# csv table')
print('{:<17}{:<15}{:<13}{:<15}{:<8}{:<10}{:<10}{:<10}{}'.format('dataset', 'sampling_type', 'tool', 'assembly_type', 'sample', 'HQ_count', 'NC_count', 'MQ_count', 'total_bin_count'))

# extract eval data
for path in RESULT_PATHS:
    # 1. dataset
    # 2. sampling type
    # 3. tool
    # (4. assembly type)
    # (5. sample)

    # extract dataset, tool, mode info from path to count to the correct dict entry
    path_list = path.split('/')

    if len(path_list) == 4:     # if coassembly append missing info
        path_list.append('NaN')
        path_list.append('all')

    dataset =       path_list[0]
    sampling_type = path_list[2]
    tool =          path_list[3]
    assembly_type = path_list[4]
    sample        = path_list[5]
    if sampling_type == 'coassembly':
        mode = 'short_co'
    else:
        mode = assembly_type + '_' + sampling_type

    # init counters per sample
    HQ_count = 0
    NC_count = 0
    MQ_count = 0
    total_bin_count = 0

    # to handle the files in that way is a Chat-GPT suggestion
    with open(path + '_check/quality_report.tsv', newline='', encoding='utf-8') as checkm2_report, \
         open(path + '/tRNA_count.txt', encoding='utf-8') as tRNA_count_txt:

        # prepare checkm2 report
        comp_cont_reader = csv.reader(checkm2_report, delimiter='\t')
        next(comp_cont_reader, None)  # skip header - also shifts start index 0 to next line(!)
        # prepare tRNA counts
        tRNA_counts = tRNA_count_txt.read().splitlines()

        for i, row in enumerate(comp_cont_reader):

            bin_ID        = row[0]  # must be of type string for use in pathnames later
            completeness  = float(row[1])
            contamination = float(row[2])
            tRNA_count    = tRNA_counts[i]  # the corresponding tRNA count


            # eval rRNA for each binning seperately, because of barrnaps output policy
            # but first handle variable fasta file namings
            fasta_add = ''
            if tool in NOT_FA_BUT_FASTA:
                fasta_add = 'sta'
            with open(path + '.' + bin_ID + '.fa' + fasta_add + '.all_rRNA_present.txt', encoding='utf-8') as rRNA_file:
                # eval if marker is set to 1 (in case all rRNA was found in this binning)
                all_rRNA_present = (rRNA_file.readline().strip() == '1')

            # count
            # increase total bin count anyway
            eval_bin_counts_dict[dataset][tool][mode]['total_bin_count'] += 1
            total_bin_count += 1
            # MQ - >50% completeness, <10% contamination
            if completeness > 0.5 and contamination < 0.1:
                eval_bin_counts_dict[dataset][tool][mode]['MQ'] += 1
                MQ_count += 1
                # append found MAG to the list of all found MAGs with every info to identify it later
                mag = (dataset, sampling_type, tool, assembly_type, sample, bin_ID)
                MAGs.append(mag)
            # NC - >90% completeness, <5% contamination
            if completeness > 0.9 and contamination < 0.05:
                eval_bin_counts_dict[dataset][tool][mode]['NC'] += 1
                NC_count += 1
                # HQ - all rRNA present?
                if all_rRNA_present:
                    eval_bin_counts_dict[dataset][tool][mode]['HQ'] += 1
                    HQ_count += 1

    print('{:<17}{:<15}{:<13}{:<15}{:<8}{:<10}{:<10}{:<10}{}'.format(dataset, sampling_type, tool, assembly_type, sample, HQ_count, NC_count, MQ_count, total_bin_count))

print('# plotting data')

for d in DATASETS:
    print('$ ' + d)
    for bq in BIN_QUALITIES:
        print('% ' + bq)
        for t in TOOLS:
            print('{:<13}{:<6}{:<6}{:<6}{:<6}{:<6}{:<6}{:<6}'.format(t, *(eval_bin_counts_dict[d][t][m][bq] for m in MODES)))
