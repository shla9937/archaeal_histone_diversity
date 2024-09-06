#!/usr/bin/env python

import sys
import statistics
from statistics import mode
#from statistics import multimode
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
from collections import Counter
from itertools import groupby
import argparse


def main():
    parser = argparse.ArgumentParser(description='Plot conservation of a protein.')
    parser.add_argument('--input', required=True, type=str,
                        help='MSA file to be plotted')
    parser.add_argument('--cutoff', required=False, type=int, default=-1,
                        help='conservation cutoff')
    parser.add_argument('--gaps', required=False, type=str, default='False',
                        help='Keep gaps? True/False')
    parser.add_argument('--cluster', required=False, type=str, default='True',
                        help='Color clustered regions? True/False')
    parser.add_argument('--sort', required=False, type=str, default='residue',
                        help='residue/mw/type')
    parser.add_argument('--color', required=False, type=str, default='orange',
                        help='Color of conserved')
    parser.add_argument('--title', required=True, type=str,
                        help='Title of plot')
    parser.add_argument('--gray', required=False, type=str, default='bbbbbb',
                        help='Which gray hexcode to use')
    parser.add_argument('--output', required=True, type=str,
                        help='Name of output file')
    parser.add_argument('--tick', required=False, type=float, default=5.5,
                        help='Size of ticks')
    parser.add_argument('--text', required=False, type=float, default=8,
                        help='Size of titles')
    args = parser.parse_args()

    input_file = str(args.input)
    conserv_cutoff = int(args.cutoff)
    keep_gaps = str(args.gaps)
    is_gradient = str(args.cluster) 
    sort_by = str(args.sort)
    color = str(args.color)
    plot_title = str(args.title)
    grey = str(args.gray)
    outfile = str(args.output)
    tick_size = args.tick
    text_size = args.text

    #open msa file from aligned mutli fasta
    msa_input_file = open(input_file, 'r', encoding="ISO-8859-1")
    sequences = import_msa(msa_input_file) 
    compiled_sequences = compile_sequences(sequences)

    #check to see if correct format, create list that is weighted for each position
    if is_msa(compiled_sequences) is True:
        compiled_positions = compile_positions(compiled_sequences)
        if sort_by == 'type':
            type_compiled_positions = sort_into_type(compiled_positions)
            hist_input, consensus = histone_type_inputer(type_compiled_positions, keep_gaps)
            conservations = get_type_conservations(type_compiled_positions, keep_gaps)
        elif sort_by == 'mw':
            mw_compiled_positions = sort_into_mw(compiled_positions)
            hist_input, consensus = histone_inputer(mw_compiled_positions, keep_gaps)
            conservations = get_conservations(mw_compiled_positions, keep_gaps)
        else:
            #make lists, remove gaps, compile consensus sequence
            hist_input, consensus = histone_inputer(compiled_positions, keep_gaps)
            #find conservation at each position
            conservations = get_conservations(compiled_positions, keep_gaps)
    else:
        print('Is not an MSA file.')
        sys.exit()
    
    consensus_str = ""
    print(consensus_str.join(consensus))
    l_curves = get_l_curve(conservations)
    if conserv_cutoff == -1:
        conserv_cutoff = int(np.mean(l_curves)+np.std(l_curves))
        if conserv_cutoff > 95:
            conserv_cutoff = 95
    #make the imput list for bar graphs
    conserv_hist_inputs = make_conservation_inputs(hist_input, conservations, conserv_cutoff)
    #make gradient inputs and plot
    if is_gradient == 'True':
        conserv_hist_inputs = make_gradient_inputs(conserv_hist_inputs)        
    elif is_gradient == 'False':
        None
    plot_histone_gram(conserv_hist_inputs, compiled_sequences, plot_title, consensus, 
                      conserv_cutoff, l_curves, color, grey, outfile, tick_size, text_size, sort_by)
    return None

def import_msa(msa):
    data = []
    sequence = None 
    seq_len = 0
    gene_count = 0

    for line in msa:
        if '>' in line:
            gene_count += 1
            if sequence is not None:
                seq_len = len(sequence)
                protein.append(sequence)
                protein.append(seq_len)
                sequence = None
                data.append(protein)
            protein = []
            protein.append(line.rstrip('\n'),)
        elif (len(line) != 0) and (sequence is not None):
            sequence += (line.rstrip('\n'))
        elif (len(line) != 0) and (sequence is None):
            sequence = (line.rstrip('\n'))
    seq_len = len(sequence)
    protein.append(sequence)
    protein.append(seq_len)
    data.append(protein)
    
    msa.close()
    return data

def compile_sequences(names_seqs_lens):
    sequences = []
    for gene in names_seqs_lens:
        sequences.append(gene[1])
    return sequences

def is_msa(sequences):
    alignment_len = -1
    sequence_idx = 0
    for sequence in sequences:
        if alignment_len == -1:
            alignment_len = len(sequence)
            sequence_idx += 1
            if len(sequences) == sequence_idx:
                return True
                break
        elif alignment_len == len(sequence):
            sequence_idx += 1
            if len(sequences) == sequence_idx:
                return True
                break
        elif alignment_len != len(sequence):
            print('Sequences are not aligned')
            return False
            break

def compile_positions(sequences):
    compiled_positions = []
    position = []
    alignment_len = len(sequences[0])
    for i in range(alignment_len):
        for sequence in sequences:
            position.append(sequence[i])
        compiled_positions.append(position)
        position = []
    return compiled_positions

def histone_inputer(compiled_positions, keep_gaps):
    hist_input = []
    consensus = []
    counter = 0 
    for i in range(len(compiled_positions)):
           #most_often = mode(compiled_positions[i])
           #most_oftens = multimode(compiled_positions[i])
        freqs = groupby(Counter(compiled_positions[i]).most_common(), lambda x:x[1])
        most_oftens = [val for val,count in next(freqs)[1]]
        most_often = most_oftens[0]

        if most_often == '-' and keep_gaps == 'False':
            continue
        else:
            for letter in compiled_positions[i]:
                if letter == most_often:
                    counter += 1
                else:
                    None
        
        for k in range(counter):
            hist_input.append(i)
        counter = 0
        consensus.append(most_often)
    return hist_input, consensus

def histone_type_inputer(compiled_positions, keep_gaps):
    hist_input = []
    consensus = []
    counter = 0
    for i in range(len(compiled_positions)):
           #most_often = mode(compiled_positions[i])
           #most_oftens = multimode(compiled_positions[i])
        freqs = groupby(Counter(compiled_positions[i]).most_common(), lambda x:x[1])
        most_oftens = [val for val,count in next(freqs)[1]]
        most_often = most_oftens[0]

        if most_often == '|' and keep_gaps == 'False':
            continue
        else:
            for letter in compiled_positions[i]:
                if letter == most_often:
                    counter += 1
                else:
                    None

        for k in range(counter):
            hist_input.append(i)
        counter = 0
        consensus.append(most_often)
    return hist_input, consensus

def get_type_conservations(compiled_positions, keep_gaps):
    conservations = []
    counter = 0
    total = len(compiled_positions[1])
    for i in range(len(compiled_positions)):
        freqs = groupby(Counter(compiled_positions[i]).most_common(), lambda x:x[1])
        most_oftens = [val for val,count in next(freqs)[1]]
        most_often = most_oftens[0]
        #most_often = mode(compiled_positions[i])
        if most_often == '|' and keep_gaps == 'False':
            None
        else:
            for letter in compiled_positions[i]:
                if letter == most_often:
                    counter += 1
                else:
                    None
        conservations.append((counter/total)*100)
        counter = 0
    return conservations

def get_conservations(compiled_positions, keep_gaps):
    conservations = []
    counter = 0
    total = len(compiled_positions[1])
    for i in range(len(compiled_positions)):
        freqs = groupby(Counter(compiled_positions[i]).most_common(), lambda x:x[1])
        most_oftens = [val for val,count in next(freqs)[1]]
        most_often = most_oftens[0]
        #most_often = mode(compiled_positions[i])
        if most_often == '-' and keep_gaps == 'False':
            None
        else:
            for letter in compiled_positions[i]:
                if letter == most_often:
                    counter += 1
                else:
                    None
        conservations.append((counter/total)*100)
        counter = 0
    return conservations

def make_conservation_inputs(hist_input, conservations, conserv_cutoff):
    conserv_idx = 0
    ungapped_idx = 0
    conserv_cutoff = conserv_cutoff
    hist_0 = []
    hist_1 = []
    conserv_hist_inputs = []
    counter = 0
    bar_info = []
    for conservation in conservations:
        if conservation >= conserv_cutoff:
            for input in hist_input:
                if input == conserv_idx:
                    counter += 1
            bar_info.append(ungapped_idx)
            bar_info.append(counter)
            counter = 0
            hist_0.append(bar_info) 
            bar_info = []
            ungapped_idx +=1   
            conserv_idx += 1
        elif conservation == 0:
            #bar_info.append(conserv_idx)
            #bar_info.append(0)
            #hist_1.append(bar_info)
            #bar_info = []
            conserv_idx += 1
            None            
        else:
            for input in hist_input:
                if input == conserv_idx:
                    counter += 1
            bar_info.append(ungapped_idx)
            bar_info.append(counter)
            counter = 0
            hist_1.append(bar_info)
            bar_info = []
            conserv_idx += 1
            ungapped_idx += 1
    conserv_hist_inputs.append(hist_0)
    conserv_hist_inputs.append(hist_1)
    return conserv_hist_inputs

def make_gradient_inputs(conserv_hist_inputs):
    cluster_5_or_over = []
    cluster_4 = []
    cluster_3 = []
    cluster_2 = []
    no_cluster = []
    below_threshold = conserv_hist_inputs[1]
    current_cluster = []
    residue_idx = 0

    for i in conserv_hist_inputs[0]:
        if i[0] == residue_idx:
            current_cluster.append(i)
            residue_idx += 1
        else:
            if len(current_cluster) >= 5:
                for j in current_cluster:
                    cluster_5_or_over.append(j)
                current_cluster = []
                current_cluster.append(i)
                residue_idx = i[0]+1
            elif len(current_cluster) == 4:
                for j in current_cluster:
                    cluster_4.append(j)
                current_cluster = []
                current_cluster.append(i)
                residue_idx = i[0]+1
            elif len(current_cluster) == 3:
                for j in current_cluster:
                    cluster_3.append(j)
                current_cluster = []
                current_cluster.append(i)
                residue_idx = i[0]+1
            elif len(current_cluster) == 2:
                for j in current_cluster:
                    cluster_2.append(j)
                current_cluster = []
                current_cluster.append(i)                
                residue_idx = i[0]+1
            elif len(current_cluster) == 1:
                for j in current_cluster:
                    no_cluster.append(j)
                current_cluster = []
                current_cluster.append(i)                
                residue_idx = i[0]+1
            elif len(current_cluster) == 0:
                current_cluster = []
                current_cluster.append(i)
                residue_idx = i[0]+1

    if len(current_cluster) >= 5:
        for j in current_cluster:
            cluster_5_or_over.append(j)
    elif len(current_cluster) == 4:
        for j in current_cluster:
            cluster_4.append(j)
    elif len(current_cluster) == 3:
        for j in current_cluster:
            cluster_3.append(j)
    elif len(current_cluster) == 2:
        for j in current_cluster:
            cluster_2.append(j)
    elif len(current_cluster) == 1:
        for j in current_cluster:
            no_cluster.append(j)

    conserv_hist_inputs = [cluster_5_or_over, cluster_4, cluster_3, cluster_2, no_cluster, below_threshold]
    return conserv_hist_inputs

def sort_into_type(compiled_positions):
    positive = ['R', 'H', 'K']
    negative = ['D', 'E']
    polar = ['S', 'T', 'N', 'Q']
    hydrophobic = ['A', 'I', 'L', 'M', 'V']
    aromatic = ['F', 'W', 'Y']
    gap = ['-']
    new_residues = []
    type_compiled_positions = []
    for sequence in compiled_positions:
        for residue in sequence:
            if residue in positive:
                new_residues.append('+')
            elif residue in negative:
                new_residues.append('-')
            elif residue in polar:
                new_residues.append('\u03B6')
            elif residue in hydrophobic:
                new_residues.append('\u03D5')
            elif residue in aromatic:
                new_residues.append('\u03A9')
            elif residue in gap:
                new_residues.append('|')
            elif residue == 'P':
                new_residues.append('P')
            elif residue == 'C':
                new_residues.append('C')
            elif residue == 'G':
                new_residues.append('G')
        type_compiled_positions.append(new_residues)
        new_residues = []

    return type_compiled_positions

def sort_into_mw(compiled_positions):
    small = ['G', 'A', 'S']
    medium_small = ['P', 'V', 'T', 'C']
    medium = ['L', 'I', 'D', 'N']
    medium_large = ['Q', 'E', 'K', 'M', 'H']
    large = ['F', 'R', 'Y', 'W']
    new_residues = []
    mw_compiled_positions = []
    for sequence in compiled_positions:
        for residue in sequence:
            if residue in small:
                new_residues.append('s')
            elif residue in medium_small:
                new_residues.append('S')
            elif residue in medium:
                new_residues.append('M')
            elif residue in medium_large:
                new_residues.append('l')
            elif residue in large:
                new_residues.append('L')
            elif residue == '-':
                new_residues.append('-')
        mw_compiled_positions.append(new_residues)
        new_residues = []

    return mw_compiled_positions

def get_l_curve(conservations):
    l_curves = []
    for i in conservations:
        if i != 0:
            l_curves.append(int(i))
    return l_curves

def get_color(color, intensity):
    oranges = ['#FBC68B', '#E69F00', '#D98D00', '#CC7B00', '#BF6A00', '#B25800']
    purples = ['#D8A2D3', '#CC79A7', '#BD6A95', '#AF5A83', '#A14A72', '#923A60']
    greens =  ['#4AC18D', '#009E73', '#008F68', '#00825D', '#007352', '#006546']
    reds =    ['#EE714B', '#D55E00', '#C15300', '#AC4800', '#973D00', '#833200']
    blues =   ['#88CAF1', '#56B4E9', '#4FA2D6', '#478FC2', '#407DAF', '#386B9B']
    yellows = ['#F4EA8C', '#F0E442', '#E6D83B', '#DBCD33', '#D1C12C', '#C7B524']
    navies =  ['#6DA2D5', '#0072B2', '#00669F', '#005B8C', '#005078', '#004564']
    blacks = ['#666666', '#555555', '#444444', '#333333', '#222222']

    if color == 'orange':
        color_shade = oranges[intensity]
    elif color == 'purple':
        color_shade = purples[intensity]
    elif color == 'green':
        color_shade = greens[intensity]
    elif color == 'red':
        color_shade = reds[intensity]
    elif color == 'blue':
        color_shade = blues[intensity]
    elif color == 'yellow':
        color_shade = yellows[intensity]
    elif color == 'navy':
        color_shade = navies[intensity]
    elif color == 'black':
        color_shade = blacks[intensity]
    return color_shade

def plot_histone_gram(conserv_hist_inputs, compiled_sequences, plot_title, consensus, 
                      conserv_cutoff, l_curves, color, grey, outfile, tick_size, text_size, sort_by):
    num_bins = len(consensus)
    n_seqs = len(compiled_sequences)
    values_high = []
    heights_high = []
    values_4 = []
    heights_4 = []
    values_3 = []
    heights_3 = []
    values_2 = []
    heights_2 = []
    values_1 = []
    heights_1 = []
    values_low = []
    heights_low = []
    xticks = [1]
    xtick_pos = [.5]
    for i in range(10, len(consensus), 10):
        xticks.append(i)
        xtick_pos.append(i-.5)
    xticks.append(None)
    xtick_pos.append(len(consensus))

    #if len(consensus) > 200:
    #    tick_size = 3
    #    fig, axs = plt.subplots(1, figsize=[num_bins*.15, 4])
    #elif len(consensus) > 150:
    #    tick_size = 5
    #    fig, axs = plt.subplots(1, figsize=[num_bins*.15, 4])
    #else:
    #    tick_size = 8
    #    fig, axs = plt.subplots(1, figsize=[num_bins*.15, 4])

    fig, axs = plt.subplots(1, figsize=[num_bins*.075, 2])

    for input in conserv_hist_inputs[0]:
        values_high.append(input[0])
        heights_high.append(input[1])
    for input in conserv_hist_inputs[-1]:
        values_low.append(input[0])
        heights_low.append(input[1])

    if len(conserv_hist_inputs) == 2:
       # axs.bar(values_high, heights_high, color=get_color(color, 5), label='>'+str(conserv_cutoff)+'% conserved')
       axs.bar(values_high, heights_high, color=get_color(color, 2), label='Highly conserved') 
       axs.bar(values_low, heights_low, color='#'+grey)
 
    elif len(conserv_hist_inputs) == 6:
        for input in conserv_hist_inputs[1]:
            values_4.append(input[0])
            heights_4.append(input[1])
        for input in conserv_hist_inputs[2]:
            values_3.append(input[0])
            heights_3.append(input[1])
        for input in conserv_hist_inputs[3]:
            values_2.append(input[0])
            heights_2.append(input[1])
        for input in conserv_hist_inputs[4]:
            values_1.append(input[0])
            heights_1.append(input[1])
        label_high = None
        label_4 = None
        label_3 = None
        label_2 = None
        label_1 = None
        label_low = None
       
        if len(values_high) > 0:
        #    label_high = '>'+str(conserv_cutoff)+'% conserved'
            label_high = 'Highly conserved'
        elif len(values_4) > 0:
        #    label_4 = '>'+str(conserv_cutoff)+'% conserved'
            label_4 = 'Highly conserved'
        elif len(values_3) > 0:
        #    label_3 = '>'+str(conserv_cutoff)+'% conserved'
            label_3 = 'Highly conserved'
        elif len(values_2) > 0:
        #    label_2 = '>'+str(conserv_cutoff)+'% conserved'
            label_2 = 'Highly conserved'
        elif len(values_1) > 0:
        #    label_1 = '>'+str(conserv_cutoff)+'% conserved'
            label_1 = 'Highly conserved'
        elif len(values_low) > 0:
        #    label_low = '<'+str(conserv_cutoff)+'% conserved'
            label_low = 'Not highly conserved' 
        axs.bar(values_high, heights_high, color=get_color(color, 4), label=label_high)
        axs.bar(values_4, heights_4, color=get_color(color, 3), label=label_4)
        axs.bar(values_3, heights_3, color=get_color(color, 2), label=label_3)
        axs.bar(values_2, heights_2, color=get_color(color, 1), label=label_2)
        axs.bar(values_1, heights_1, color=get_color(color, 0), label=label_1)
        axs.bar(values_low, heights_low, color='#'+grey, label=label_low)    

    axs.legend(loc='lower right', bbox_to_anchor=(1, 1), fontsize=tick_size, markerscale=0.5)
    plt.title(plot_title+'\n (n='+str(n_seqs)+')', fontsize=text_size)
    plt.xlabel('Residue'+' (len='+str(len(consensus))+')', labelpad=3*tick_size, fontsize=text_size)
    plt.ylabel('Conservation', fontsize=text_size)
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    percent_max = len(compiled_sequences)
    axs.set_ylim(0,percent_max)
    axs.yaxis.set_major_formatter(PercentFormatter(xmax=percent_max))
    plt.yticks([0, n_seqs*.25, n_seqs*.5, n_seqs*.75, n_seqs], ['0%', '25%', '50%', '75%', '100%'], fontsize=tick_size)
    #if sort_by == 'type':
    #    plt.xticks(range(len(consensus)), consensus, size=tick_size, rotation='vertical')
    #else:
    plt.xticks(range(len(consensus)), consensus, size=tick_size)
    plt.tick_params(axis='x', pad=.5*tick_size, top=False, bottom=False)
    plt.plot()  
    
    axs_2 = plt.twiny()
    axs_2.spines['right'].set_visible(False)
    axs_2.spines['top'].set_visible(False)
    #if sort_by == 'type':
    #    plt.tick_params(axis='x', labelbottom=True, labeltop=False, top=False, bottom=False, pad=4*tick_size)
    #else:
    plt.tick_params(axis='x', labelbottom=True, labeltop=False, top=False, bottom=False, pad=2*tick_size)
    
    plt.xticks(xtick_pos, xticks, size=tick_size)
    axs.margins(y=0, x=0)
    axs_2.margins(y=0, x=0)
    plt.tight_layout()
    plt.savefig(outfile, format='pdf')
    plt.plot()    

    #plt.figure(2)
    #plt.axvline(x=np.mean(l_curves)-np.std(l_curves), ls = "--", color='#aaaaaa', alpha=0.7)
    #plt.axvline(x=np.mean(l_curves)+np.std(l_curves), ls = "--", color='#aaaaaa', alpha=0.7)
    #mu = np.mean(l_curves)
    #sigma = np.std(l_curves)
    #n, bins, patches = plt.hist(l_curves, density=1)
    #y = ((1 / (np.sqrt(2 * np.pi) * sigma)) *
    #    np.exp(-0.5 * (1 / sigma * (bins - mu))**2))
    #plt.plot(bins, y, '--')
    plt.show()


if __name__ == '__main__':
    main()

