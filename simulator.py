#! /usr/bin/python
# -*- coding: utf-8 -*-

import random
import argparse
import re
from lib import Fasta, write_fastas_to_a_file

FASTA_PATTERN = re.compile('>(\d+)\s+([^>]+)', re.DOTALL)


def read_fasta_file(path):
    with open(path) as f:
        return [Fasta(id_, ''.join(seq.split())) for id_, seq in re.findall(FASTA_PATTERN, f.read())]


def create_random_fragments_for_set_of_seqs(fasta_sequences, length, num_of_fragments, substitution_rate,
                                            possible_base_symbols):
    fragments = []
    length_of_start_fragment_for_a_sequence = {sequence: len(sequence) - length for sequence in fasta_sequences if
                                               len(sequence) - length >= 0}
    intervals_for_sequences = {}
    temp_sum = 0
    for seq, length_of_start_fragment in length_of_start_fragment_for_a_sequence.iteritems():
        intervals_for_sequences[(temp_sum, temp_sum + length_of_start_fragment + 1)] = seq
        temp_sum += length_of_start_fragment + 1
    while len(fragments) < num_of_fragments:
        x = random.randint(0, temp_sum - 1)
        interval = get_interval(x, intervals_for_sequences.keys())
        seq = intervals_for_sequences[interval].seq
        description = intervals_for_sequences[interval].description
        fragments.append(Fasta(description, create_random_fragment(seq, length, substitution_rate,possible_base_symbols)))
    return fragments


def get_interval(number, intervals):
    for interval in intervals:
        if interval[0] <= number < interval[1]:
            return interval


def substitute_base(sequence, substitution_rate, possible_values):
    if random.random() < substitution_rate:
        position = random.randint(0, len(sequence) - 1)
        possible_bases = [val for val in possible_values if val != sequence[position]]
        base = random.choice(possible_bases)
        return sequence[:position] + base + sequence[position + 1:]
    return sequence


def create_random_fragment(sequence, length, substitution_rate, possible_base_symbols):
    assert len(sequence) >= length
    start_max = len(sequence) - length
    start = random.randint(0, start_max)
    return substitute_base(sequence[start:start + length], substitution_rate, possible_base_symbols)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')  # todo dodać opis
    parser.add_argument('sequence', type=str, help='sequence file in FASTA format')
    parser.add_argument('length', type=int, help='length of the obtained fragments')
    parser.add_argument('num_of_fragments', type=int, help='number of fragments to be created')
    parser.add_argument('substitution_rate', type=float, help='probability of substitution of a base')
    parser.add_argument('outfile', type=str, help='path to the outfile')
    parser.add_argument('--possible_bases', type=str, default='ACGT')
    parser.add_argument('--seed', type=int, default=77)
    args = parser.parse_args()
    assert args.length > 0
    assert args.num_of_fragments > 0
    assert 0 <= args.substitution_rate <= 1
    random.seed(args.seed)
    fasta_seqs = read_fasta_file(args.sequence)
    print "Read %d sequences" % len(fasta_seqs)
    fasta_fragments = create_random_fragments_for_set_of_seqs(fasta_seqs, args.length, args.num_of_fragments, args.substitution_rate, args.possible_bases)
    write_fastas_to_a_file(fasta_fragments, args.outfile)
    print "Fragments written to", args.outfile
