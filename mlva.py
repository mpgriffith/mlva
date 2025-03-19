#!/usr/bin/env python

import argparse
import sys
import os
import re
from io import StringIO
import csv
import regex

from Bio import SeqIO, SearchIO
from Bio.Blast.Applications import NcbiblastnCommandline as blastn

DEFAULT_OUTFMT = '"6 std qlen slen stitle qcovs qframe sframe nident"'
DEFAULT_FIELDS = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
        'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
extra_fields = ['qlen', 'slen', 'stitle', 'qcovs', 'qframe', 'sframe', 'nident']


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fastas', nargs='+')
    parser.add_argument('-v', '--primer_fasta')
    parser.add_argument('-o', '--output', default=None)
    parser.add_argument('-i', '--verbose_output', default=None)
    parser.add_argument('-c', '--count', default=False, action='store_true')

    return parser.parse_args()

def find_repeats(seq, repeat_seq):
    # print(seq, repeat_seq)
    er_num = 1 #str(int(len(repeat_seq) /2))
    # patt = '((?<=%s)%s)' % (repeat_seq, repeat_seq)
    pattern = '(%s){e<=%s}' % (repeat_seq, str(1))
    #pattern = '(%s)' % repeat_seq
    # print(pattern)
    repeats = regex.findall(pattern, str(seq.seq), flags=regex.IGNORECASE)
    # print(pattern, len(repeats))
    if len(repeats) == 0:
        repeats = regex.findall(pattern, str(seq.reverse_complement().seq), flags=regex.IGNORECASE)
    # print(len(repeats))
    return len(repeats) 

def find_all_vntrs(fasta_file, primer_file, vntrs, save_files=False):
    blast_cmd = blastn(task='blastn-short', query=primer_file, subject=fasta_file, outfmt=DEFAULT_OUTFMT)
    bc_out, bc_err = blast_cmd()
    bf = StringIO(bc_out)
    blast_results = SearchIO.parse(bf, 'blast-tab', fields=DEFAULT_FIELDS+extra_fields)
    best_primer_hits = {}
    for qr in blast_results:
        # print(qr)
        # qr = qr.sort(key=lambda h: (-h[0].evalue, h[0].ident_pct))
        # print(qr)
        best_primer_hits[qr.id] = qr if qr and len(qr) > 0 and len(qr[0]) > 0 else None

    data = {v: None for v in vntrs}
    fasta_seqs = SeqIO.index(fasta_file, 'fasta')
    n = fasta_file.split('/')[-1].split('.')[0]
    if save_files:
        bf = open(os.path.join(save_files, '%s_primer_blast.tsv' % n), 'w')
        bf.write('\t'.join(DEFAULT_FIELDS+extra_fields) + '\n')
        SearchIO.write([qr.hsp_filter(lambda hit: hit == qr[0][0]) for qr in best_primer_hits.values()],
         bf, 'blast-tab', fields=DEFAULT_FIELDS+extra_fields)
        bf.close()
    amp_seqs = []
    lens = {}
    for vntr, repeat in vntrs.items():
        forward_hit = best_primer_hits[vntr + 'F'][0][0] if vntr + 'F' in best_primer_hits else None
        reverse_hit = best_primer_hits[vntr + 'R'][0][0] if vntr + 'R' in best_primer_hits else None
        if not forward_hit or not reverse_hit:
            data[vntr] = 'NP'
            continue
        forward_contig = forward_hit.hit_id
        reverse_contig = reverse_hit.hit_id
        if reverse_contig != forward_contig:
            fs, fe = tuple(sorted([forward_hit.hit_start, forward_hit.hit_end]))
            forward_seq = fasta_seqs[forward_contig][fs:] if forward_hit.hit_strand > 0 else fasta_seqs[forward_contig][:fe]
            rs, re = tuple(sorted([reverse_hit.hit_start, reverse_hit.hit_end]))
            reverse_seq = fasta_seqs[reverse_contig][rs:] if reverse_hit.hit_strand > 0 else fasta_seqs[reverse_contig][:rs]
            f_repeats = find_repeats(forward_seq, vntrs[vntr])
            r_repeats = find_repeats(reverse_seq, vntrs[vntr])
            data[vntr] = '>=' + str(max([f_repeats, r_repeats]))
            continue
        forward_start, forward_end = forward_hit.hit_start, forward_hit.hit_end
        reverse_start, reverse_end = reverse_hit.hit_start, reverse_hit.hit_end
        min_pos = min(forward_start, forward_end, reverse_start, reverse_end)
        max_pos = max(forward_start, forward_end, reverse_start, reverse_end)
        
        if reverse_start < forward_start:
            # print(reverse_start, reverse_end, forward_start, forward_end)
            seq = fasta_seqs[forward_contig][min_pos: max_pos]
        else:
            seq = fasta_seqs[forward_contig][min_pos:max_pos]
        # print(vntrs[vntr], seq.seq[10:])
        data[vntr] = find_repeats(seq, vntrs[vntr])
        if save_files:
            seq.id = '%s_%s_%s__%s' % (forward_contig, forward_start, reverse_end, vntr)
            lens[vntr] = len(seq)
            amp_seqs.append(seq)
        # print(vntr, data[vntr], len(seq))
    if save_files:
        SeqIO.write(amp_seqs, os.path.join(save_files, '%s_seqs.fasta' % n), 'fasta')
    return data, lens

def count_motif(seqs, repeat_seq):
    # TG[AC][CG] = len(s) - s.count('[') * 3
    max_num_repeats = 0
    num_errors_per_repeat = 1#max(1, int((len(repeat_seq) - (repeat_seq.count('[') * 3)) / 2 ) - 3)
    for s in seqs:
        i = 0
        found = True
        while found:
            motif_pattern = r'(%s){e<=%s}' % (repeat_seq * (i+1), str(num_errors_per_repeat * (i+1)))
            # print(motif_pattern)
            forward_hits = regex.findall(motif_pattern, str(s.seq), flags=regex.IGNORECASE)
            reverse_hits = regex.findall(motif_pattern, str(s.seq.reverse_complement()), flags=regex.IGNORECASE)
            found = True if len(forward_hits) > 0 or len(reverse_hits) > 0 else False
            i += 1
        if i > max_num_repeats: max_num_repeats = i
    return max_num_repeats


def count_all_repeats(fasta_file, vntrs):
    seqs = list(SeqIO.parse(fasta_file, 'fasta'))
    data = {v: None for v in vntrs}
    for v, r in vntrs.items():
        data[v] = count_motif(seqs, r)
    return data, None


def find_repeats_by_length(fasta_file, primer_file, vntrs, save_files=False):
    blast_cmd = blastn(task='blastn-short', query=primer_file, subject=fasta_file, outfmt=DEFAULT_OUTFMT)
    bc_out, bc_err = blast_cmd()
    bf = StringIO(bc_out)
    blast_results = SearchIO.parse(bf, 'blast-tab', fields=DEFAULT_FIELDS+extra_fields)
    best_primer_hits = {}
    for qr in blast_results:
        # print(qr)
        # qr = qr.sort(key=lambda h: (-h[0].evalue, h[0].ident_pct))
        # print(qr)
        best_primer_hits[qr.id] = qr if qr and len(qr) > 0 and len(qr[0]) > 0 else None
    vntr_lens = {s.id[:-1]: int(s.description.split()[-1]) for s in SeqIO.parse(primer_file, 'fasta')}
    data = {v: None for v in vntrs}
    fasta_seqs = SeqIO.index(fasta_file, 'fasta')
    n = fasta_file.split('/')[-1].split('.')[0]
    if save_files:
        bf = open(os.path.join(save_files, '%s_primer_blast.tsv' % n), 'w')
        bf.write('\t'.join(DEFAULT_FIELDS+extra_fields) + '\n')
        SearchIO.write([qr.hsp_filter(lambda hit: hit == qr[0][0]) for qr in best_primer_hits.values()],
         bf, 'blast-tab', fields=DEFAULT_FIELDS+extra_fields)
        bf.close()
        with open(os.path.join(save_files, '%s_raw_primer_blast.tsv' % n), 'w') as raw_b:
             raw_b.write(bc_out)

    amp_seqs = []
    lens = {}
    for vntr, repeat in vntrs.items():
        forward_hit = best_primer_hits[vntr + 'F'][0][0] if vntr + 'F' in best_primer_hits else None
        reverse_hit = best_primer_hits[vntr + 'R'][0][0] if vntr + 'R' in best_primer_hits else None
        if not forward_hit or not reverse_hit:
            data[vntr] = 'NP'
            continue
        forward_contig = forward_hit.hit_id
        reverse_contig = reverse_hit.hit_id
        if reverse_contig != forward_contig:
            data[vntr] = 'DC'
            continue
        forward_start, forward_end = forward_hit.hit_start, forward_hit.hit_end
        reverse_start, reverse_end = reverse_hit.hit_start, reverse_hit.hit_end
        min_pos = min(forward_start, forward_end, reverse_start, reverse_end)
        max_pos = max(forward_start, forward_end, reverse_start, reverse_end)
        
        if reverse_start < forward_start:
            # print(reverse_start, reverse_end, forward_start, forward_end)
            seq = fasta_seqs[forward_contig][min_pos: max_pos]
        else:
            seq = fasta_seqs[forward_contig][min_pos:max_pos]
        # print(vntrs[vntr], seq.seq[10:])
        repeat_len = len(vntrs[vntr]) - vntrs[vntr].count('[') * 3
        data[vntr] = round((len(seq) - vntr_lens[vntr])/repeat_len)
        if save_files:
            seq.id = '%s_%s_%s__%s' % (forward_contig, forward_start, reverse_end, vntr)
            lens[vntr] = len(seq)
            amp_seqs.append(seq)
        # print(vntr, data[vntr], len(seq))
    if save_files:
        SeqIO.write(amp_seqs, os.path.join(save_files, '%s_seqs.fasta' % n), 'fasta')
    return data, lens

def main():
    
    args = parse_args()
    fastas = args.fastas
    verbose_output = args.verbose_output
    verbose = verbose_output if verbose_output else False
    primer_file = args.primer_fasta
    if len(fastas) == 1 and os.path.splitext(fastas[0])[1][1:] not in ['fasta', 'faa', 'fa', 'fsa', 'fna']:
        fasta_files = []
        with open(fastas[0], 'r') as fasta_list_file:
            for line in fasta_list_file:
                fasta_files.append(line.strip())
        fastas = fasta_files
    vntrs = {}
    primers = SeqIO.parse(primer_file, 'fasta')
    for p in primers:
        seq_id = p.id
        v = seq_id[0:-1]
        vntrs[v] = p.description.split()[1].strip()
    vntr_data = {}
    lens = {}
    # print(fastas)
    for f in fastas:
        n = f.split('/')[-1].rsplit('.')[0]
        if args.count:
            vntr_data[n], lens[n] = find_all_vntrs(f, primer_file, vntrs, save_files=verbose)
        else:
            vntr_data[n], lens[n] = find_repeats_by_length(f, primer_file, vntrs, save_files=verbose)
        # vntr_data[n], _ = find_all_vntrs(f, primer_file, vntrs)
     # print(vntrs)
    output = open(args.output, 'w') if args.output else sys.stdout
    writer = csv.DictWriter(output, ['Isolate'] + list(vntrs.keys()))
    writer.writeheader()
    for n, d in vntr_data.items():
        d['Isolate'] = n
        # print(d)
        writer.writerow(d)
    if verbose:
        with open(os.path.join(verbose, 'vntr_lengths.csv'), 'w') as f:
            writer = csv.DictWriter(f, ['Isolate'] + list(vntrs.keys()))
            writer.writeheader()
            for n, l in lens.items():
                l['Isolate'] = n
                writer.writerow(l)
    

if __name__ == "__main__":
    main()
