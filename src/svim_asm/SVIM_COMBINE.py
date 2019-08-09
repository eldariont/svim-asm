import os
import logging
import re

from collections import defaultdict
from math import pow, sqrt
import time
from statistics import mean, stdev
from edlib import align
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster

from svim_asm.SVCandidate import CandidateInversion, CandidateDuplicationTandem, CandidateDuplicationInterspersed, CandidateDeletion, CandidateInsertion, CandidateBreakend

def add_reverse_of_breakends(breakend_candidates):
    final_candidates = []
    for candidate in breakend_candidates:
        final_candidates.append(candidate)
        new_source_direction = 'fwd' if candidate.dest_direction == 'rev' else 'rev'
        new_dest_direction = 'fwd' if candidate.source_direction == 'rev' else 'rev'
        final_candidates.append(CandidateBreakend(candidate.dest_contig, 
                                                  candidate.dest_start, 
                                                  new_source_direction, 
                                                  candidate.source_contig, 
                                                  candidate.source_start, 
                                                  new_dest_direction, 
                                                  candidate.reads, 
                                                  candidate.genotype))
    return final_candidates


def form_partitions(sv_candidates_with_haplotype, max_distance):
    """Form partitions of signatures using mean distance."""
    sorted_candidates_with_haplotype = sorted(sv_candidates_with_haplotype, key=lambda evi: evi[1].get_key())
    partitions = []
    current_partition = []
    for haplotype, candidate in sorted_candidates_with_haplotype:
        if len(current_partition) > 0:
            candidate_key = candidate.get_key()
            last_candidate_key = current_partition[-1][1].get_key()
            if last_candidate_key[0] != candidate_key[0] or \
               last_candidate_key[1] != candidate_key[1] or \
               abs(last_candidate_key[2] - candidate_key[2]) > max_distance:
                partitions.append(current_partition[:])
                current_partition = []
        current_partition.append((haplotype, candidate))
    if len(current_partition) > 0:
        partitions.append(current_partition[:])
    return partitions


def compute_distance(candidate_with_haplotype1, candidate_with_haplotype2, reference):
    haplotype1, candidate1 = candidate_with_haplotype1
    haplotype2, candidate2 = candidate_with_haplotype2
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    if haplotype1 == haplotype2:
        return 1000000000

    if candidate1.type == "DEL":
        region_chr = candidate1.source_contig
        region_start = min(candidate1.source_start, candidate2.source_start) - 100
        region_end = max(candidate1.source_end, candidate2.source_end) + 100
        haplotype1 = reference.fetch(region_chr, region_start, candidate1.source_start).upper() + reference.fetch(region_chr, candidate1.source_end, region_end).upper()
        haplotype2 = reference.fetch(region_chr, region_start, candidate2.source_start).upper() + reference.fetch(region_chr, candidate2.source_end, region_end).upper()
        editDistance = align(haplotype1, haplotype2)["editDistance"]
    elif candidate1.type == "INV":
        region_chr = candidate1.source_contig
        region_start = min(candidate1.source_start, candidate2.source_start) - 100
        region_end = max(candidate1.source_end, candidate2.source_end) + 100
        inverted_seq1 = "".join(complement.get(base.upper(), base.upper()) for base in reversed(reference.fetch(region_chr, candidate1.source_start, candidate1.source_end).upper()))
        haplotype1 = reference.fetch(region_chr, region_start, candidate1.source_start).upper() + \
                     inverted_seq1 + \
                     reference.fetch(region_chr, candidate1.source_end, region_end).upper()
        inverted_seq2 = "".join(complement.get(base.upper(), base.upper()) for base in reversed(reference.fetch(region_chr, candidate2.source_start, candidate2.source_end).upper()))
        haplotype2 = reference.fetch(region_chr, region_start, candidate2.source_start).upper() + \
                     inverted_seq2 + \
                     reference.fetch(region_chr, candidate2.source_end, region_end).upper()
        editDistance = align(haplotype1, haplotype2)["editDistance"]
    elif candidate1.type == "INS":
        region_chr = candidate1.dest_contig
        region_start = min(candidate1.dest_start, candidate2.dest_start) - 100
        region_end = max(candidate1.dest_start, candidate2.dest_start) + 100
        haplotype1 = reference.fetch(region_chr, region_start, candidate1.dest_start).upper() + \
                     candidate1.sequence + \
                     reference.fetch(region_chr, candidate1.dest_start, region_end).upper()
        haplotype2 = reference.fetch(region_chr, region_start, candidate2.dest_start).upper() + \
                     candidate2.sequence + \
                     reference.fetch(region_chr, candidate2.dest_start, region_end).upper()
        editDistance = align(haplotype1, haplotype2)["editDistance"]
    elif candidate1.type == "DUP_TAN":
        region_chr = candidate1.source_contig
        region_start = min(candidate1.source_start, candidate2.source_start) - 100
        region_end = max(candidate1.source_end, candidate2.source_end) + 100
        haplotype1 = reference.fetch(region_chr, region_start, candidate1.source_start).upper() + \
                     reference.fetch(region_chr, candidate1.source_start, candidate1.source_end).upper() * candidate1.copies + \
                     reference.fetch(region_chr, candidate1.source_end, region_end).upper()
        haplotype2 = reference.fetch(region_chr, region_start, candidate2.source_start).upper() + \
                     reference.fetch(region_chr, candidate2.source_start, candidate2.source_end).upper() * candidate2.copies + \
                     reference.fetch(region_chr, candidate2.source_end, region_end).upper()
        editDistance = align(haplotype1, haplotype2)["editDistance"]
    elif candidate1.type == "DUP_INT":
        region_chr = candidate1.dest_contig
        region_start = min(candidate1.dest_start, candidate2.dest_start) - 100
        region_end = max(candidate1.dest_start, candidate2.dest_start) + 100
        haplotype1 = reference.fetch(region_chr, region_start, candidate1.dest_start).upper() + \
                     reference.fetch(candidate1.source_contig, candidate1.source_start, candidate1.source_end).upper() + \
                     reference.fetch(region_chr, candidate1.dest_start, region_end).upper()
        haplotype2 = reference.fetch(region_chr, region_start, candidate2.dest_start).upper() + \
                     reference.fetch(candidate2.source_contig, candidate2.source_start, candidate2.source_end).upper() + \
                     reference.fetch(region_chr, candidate2.dest_start, region_end).upper()
        editDistance = align(haplotype1, haplotype2)["editDistance"]
    # elif candidate1.type == "BND":
    #     region_start = min(candidate1.dest_start, candidate2.dest_start) - 100
    #     region_end = max(candidate1.dest_start, candidate2.dest_start) + 100
    #     haplotype1 = reference.fetch(region_chr, region_start, candidate1.dest_start).upper() + \
    #                  reference.fetch(candidate1.source_contig, candidate1.source_start, candidate1.source_end).upper() + \
    #                  reference.fetch(region_chr, candidate1.dest_start, region_end).upper()
    #     haplotype2 = reference.fetch(region_chr, region_start, candidate2.dest_start).upper() + \
    #                  reference.fetch(candidate2.source_contig, candidate2.source_start, candidate2.source_end).upper() + \
    #                  reference.fetch(region_chr, candidate2.dest_start, region_end).upper()
    #     editDistance = align(haplotype1, haplotype2)["editDistance"]

    return editDistance


def pair_haplotypes(partitions, reference):
    clusters_final = []
    for partition in partitions:
        if len(partition) < 2:
            new_clusters = [partition]
        else:
            distances = []
            for i in range(len(partition)-1):
                for j in range(i+1, len(partition)):
                    distances.append(compute_distance(partition[i], partition[j], reference))
            Z = linkage(np.array(distances), method = "complete")
            cluster_indices = list(fcluster(Z, 10, criterion='distance'))
            new_clusters = [[] for i in range(max(cluster_indices))]
            for candidate_index, cluster_index in enumerate(cluster_indices):
                new_clusters[cluster_index-1].append(partition[candidate_index])
        clusters_final.extend(new_clusters)
    return clusters_final


def pair_candidates(sv_candidates1, sv_candidates2, reference):
    deletion_candidates1 = [(1, cand) for cand in sv_candidates1 if cand.type == "DEL"]
    insertion_candidates1 = [(1, cand) for cand in sv_candidates1 if cand.type == "INS"]
    inversion_candidates1 = [(1, cand) for cand in sv_candidates1 if cand.type == "INV"]
    tandem_duplication_candidates1 = [(1, cand) for cand in sv_candidates1 if cand.type == "DUP_TAN"]
    breakend_candidates1 = [(1, cand) for cand in sv_candidates1 if cand.type == "BND"]
    interspersed_duplication_candidates1 = [(1, cand) for cand in sv_candidates1 if cand.type == "DUP_INT"]

    deletion_candidates2 = [(2, cand) for cand in sv_candidates2 if cand.type == "DEL"]
    insertion_candidates2 = [(2, cand) for cand in sv_candidates2 if cand.type == "INS"]
    inversion_candidates2 = [(2, cand) for cand in sv_candidates2 if cand.type == "INV"]
    tandem_duplication_candidates2 = [(2, cand) for cand in sv_candidates2 if cand.type == "DUP_TAN"]
    breakend_candidates2 = [(2, cand) for cand in sv_candidates2 if cand.type == "BND"]
    interspersed_duplication_candidates2 = [(2, cand) for cand in sv_candidates2 if cand.type == "DUP_INT"]

    paired_candidates = []
    #DELETIONS
    logging.info("Pairing {0} deletions...".format(len(deletion_candidates1) + len(deletion_candidates2)))
    partitions = form_partitions(deletion_candidates1 + deletion_candidates2, 10000)
    clusters = pair_haplotypes(partitions, reference)
    for cluster in clusters:
        if len(cluster) == 1:
            candidate = cluster[0][1]
            genotype = "1/0" if cluster[0][0] == 1 else "0/1"
            paired_candidates.append(CandidateDeletion(candidate.source_contig, 
                                                       candidate.source_start, 
                                                       candidate.source_end, 
                                                       candidate.reads, 
                                                       genotype))
        elif len(cluster) == 2:
            candidate = cluster[0][1]
            reads = cluster[0][1].reads + cluster[1][1].reads
            genotype = "1/1"
            paired_candidates.append(CandidateDeletion(candidate.source_contig, 
                                                       candidate.source_start, 
                                                       candidate.source_end, 
                                                       reads, 
                                                       genotype))
        else:
            logging.error("Cluster size should be either 1 or 2 but is " + str(len(cluster)))
    
    #INVERSIONS
    logging.info("Pairing {0} inversions...".format(len(inversion_candidates1) + len(inversion_candidates2)))
    partitions = form_partitions(inversion_candidates1 + inversion_candidates2, 10000)
    clusters = pair_haplotypes(partitions, reference)
    for cluster in clusters:
        if len(cluster) == 1:
            candidate = cluster[0][1]
            genotype = "1/0" if cluster[0][0] == 1 else "0/1"
            paired_candidates.append(CandidateInversion(candidate.source_contig, 
                                                        candidate.source_start, 
                                                        candidate.source_end, 
                                                        candidate.reads, 
                                                        candidate.complete, 
                                                        genotype))
        elif len(cluster) == 2:
            candidate = cluster[0][1]
            complete = cluster[0][1].complete or cluster[1][1].complete
            reads = cluster[0][1].reads + cluster[1][1].reads
            genotype = "1/1"
            paired_candidates.append(CandidateInversion(candidate.source_contig, 
                                                        candidate.source_start, 
                                                        candidate.source_end, 
                                                        reads, 
                                                        complete, 
                                                        genotype))
        else:
            logging.error("Cluster size should be either 1 or 2 but is " + str(len(cluster)))

    #INSERTIONS
    logging.info("Pairing {0} insertions...".format(len(insertion_candidates1) + len(insertion_candidates2)))
    partitions = form_partitions(insertion_candidates1 + insertion_candidates2, 10000)
    clusters = pair_haplotypes(partitions, reference)
    for cluster in clusters:
        if len(cluster) == 1:
            candidate = cluster[0][1]
            genotype = "1/0" if cluster[0][0] == 1 else "0/1"
            paired_candidates.append(CandidateInsertion(candidate.dest_contig, 
                                                        candidate.dest_start, 
                                                        candidate.dest_end, 
                                                        candidate.reads, 
                                                        candidate.sequence, 
                                                        genotype))
        elif len(cluster) == 2:
            candidate = cluster[0][1]
            reads = cluster[0][1].reads + cluster[1][1].reads
            genotype = "1/1"
            paired_candidates.append(CandidateInsertion(candidate.dest_contig, 
                                                        candidate.dest_start, 
                                                        candidate.dest_end, 
                                                        reads, 
                                                        candidate.sequence, 
                                                        genotype))
        else:
            logging.error("Cluster size should be either 1 or 2 but is " + str(len(cluster)))

    #TANDEM DUPLICATIONS
    logging.info("Pairing {0} tandem duplications...".format(len(tandem_duplication_candidates1) + len(tandem_duplication_candidates2)))
    partitions = form_partitions(tandem_duplication_candidates1 + tandem_duplication_candidates2, 10000)
    clusters = pair_haplotypes(partitions, reference)
    for cluster in clusters:
        if len(cluster) == 1:
            candidate = cluster[0][1]
            genotype = "1/0" if cluster[0][0] == 1 else "0/1"
            paired_candidates.append(CandidateDuplicationTandem(candidate.source_contig, 
                                                                candidate.source_start, 
                                                                candidate.source_end, 
                                                                candidate.copies, 
                                                                candidate.fully_covered, 
                                                                candidate.reads,
                                                                genotype))
        elif len(cluster) == 2:
            candidate = cluster[0][1]
            fully_covered = cluster[0][1].fully_covered or cluster[1][1].fully_covered
            reads = cluster[0][1].reads + cluster[1][1].reads
            genotype = "1/1"
            paired_candidates.append(CandidateDuplicationTandem(candidate.source_contig, 
                                                                candidate.source_start, 
                                                                candidate.source_end, 
                                                                round(mean([cluster[0][1].copies, cluster[1][1].copies])),
                                                                fully_covered,
                                                                reads, 
                                                                genotype))
        else:
            logging.error("Cluster size should be either 1 or 2 but is " + str(len(cluster)))

    #INTERSPERSED DUPLICATIONS
    logging.info("Pairing {0} interspersed duplications...".format(len(interspersed_duplication_candidates1) + len(interspersed_duplication_candidates2)))
    partitions = form_partitions(interspersed_duplication_candidates1 + interspersed_duplication_candidates2, 10000)
    clusters = pair_haplotypes(partitions, reference)
    for cluster in clusters:
        if len(cluster) == 1:
            candidate = cluster[0][1]
            genotype = "1/0" if cluster[0][0] == 1 else "0/1"
            paired_candidates.append(CandidateDuplicationInterspersed(candidate.source_contig, 
                                                                      candidate.source_start, 
                                                                      candidate.source_end, 
                                                                      candidate.dest_contig, 
                                                                      candidate.dest_start, 
                                                                      candidate.dest_end,
                                                                      candidate.reads,
                                                                      candidate.cutpaste,
                                                                      genotype))
        elif len(cluster) == 2:
            candidate = cluster[0][1]
            fully_covered = cluster[0][1].fully_covered or cluster[1][1].fully_covered
            reads = cluster[0][1].reads + cluster[1][1].reads
            cutpaste = cluster[0][1].cutpaste or cluster[1][1].cutpaste
            genotype = "1/1"
            paired_candidates.append(CandidateDuplicationInterspersed(candidate.source_contig, 
                                                                      candidate.source_start, 
                                                                      candidate.source_end, 
                                                                      candidate.dest_contig, 
                                                                      candidate.dest_start, 
                                                                      candidate.dest_end,
                                                                      reads,
                                                                      cutpaste,
                                                                      genotype))
        else:
            logging.error("Cluster size should be either 1 or 2 but is " + str(len(cluster)))

    # #BREAKENDS
    # partitions = form_partitions(breakend_candidates1 + breakend_candidates2, 10000)
    # clusters = pair_haplotypes(partitions, reference)
    # for cluster in clusters:
    #     if len(cluster) == 1:
    #         candidate = cluster[0][1]
    #         genotype = "1/0" if cluster[0][0] == 1 else "0/1"
    #         paired_candidates.append(CandidateDuplicationInterspersed(candidate.source_contig, 
    #                                                                   candidate.source_start, 
    #                                                                   candidate.source_end, 
    #                                                                   candidate.dest_contig, 
    #                                                                   candidate.dest_start, 
    #                                                                   candidate.dest_end,
    #                                                                   candidate.reads,
    #                                                                   candidate.cutpaste,
    #                                                                   genotype))
    #     elif len(cluster) == 2:
    #         candidate = cluster[0][1]
    #         fully_covered = cluster[0][1].fully_covered or cluster[1][1].fully_covered
    #         reads = cluster[0][1].reads + cluster[1][1].reads
    #         cutpaste = cluster[0][1].cutpaste or cluster[1][1].cutpaste
    #         genotype = "1/1"
    #         paired_candidates.append(CandidateDuplicationInterspersed(candidate.source_contig, 
    #                                                                   candidate.source_start, 
    #                                                                   candidate.source_end, 
    #                                                                   candidate.dest_contig, 
    #                                                                   candidate.dest_start, 
    #                                                                   candidate.dest_end,
    #                                                                   reads,
    #                                                                   cutpaste,
    #                                                                   genotype))
    #     else:
    #         logging.error("Cluster size should be either 1 or 2 but is " + str(len(cluster)))
    return paired_candidates


def sorted_nicely(vcf_entries):
    """ Sort the given vcf entries (in the form ((contig, start, end), vcf_string, sv_type)) in the way that humans expect.
        e.g. chr10 comes after chr2
        Algorithm adapted from https://blog.codinghorror.com/sorting-for-humans-natural-sort-order/"""
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    tuple_key = lambda entry: ( alphanum_key(str(entry[0][0])), entry[0][1], entry[0][2] )
    return sorted(vcf_entries, key = tuple_key)


def write_final_vcf(int_duplication_candidates,
                    inversion_candidates, 
                    tandem_duplication_candidates, 
                    deletion_candidates, 
                    insertion_candidates, 
                    breakend_candidates, 
                    version, 
                    contig_names, 
                    contig_lengths,
                    types_to_output,
                    reference,
                    options):
    vcf_output = open(options.working_dir + '/variants.vcf', 'w')

    # Write header lines
    print("##fileformat=VCFv4.2", file=vcf_output)
    print("##fileDate={0}".format(time.strftime("%Y-%m-%d|%I:%M:%S%p|%Z|%z")), file=vcf_output)
    print("##source=SVIM-asm-v{0}".format(version), file=vcf_output)
    for contig_name, contig_length in zip(contig_names, contig_lengths):
        print("##contig=<ID={0},length={1}>".format(contig_name, contig_length), file=vcf_output)
    if "DEL" in types_to_output:
        print("##ALT=<ID=DEL,Description=\"Deletion\">", file=vcf_output)
    if "INV" in types_to_output:
        print("##ALT=<ID=INV,Description=\"Inversion\">", file=vcf_output)
    if not options.duplications_as_insertions and ("DUP_TAN" in types_to_output or "DUP_INT" in types_to_output):
        print("##ALT=<ID=DUP,Description=\"Duplication\">", file=vcf_output)
    if not options.duplications_as_insertions and "DUP_TAN" in types_to_output:
        print("##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">", file=vcf_output)
    if not options.duplications_as_insertions and "DUP_INT" in types_to_output:
        print("##ALT=<ID=DUP:INT,Description=\"Interspersed Duplication\">", file=vcf_output)
    if "INS" in types_to_output:
        print("##ALT=<ID=INS,Description=\"Insertion\">", file=vcf_output)
    if "BND" in types_to_output:
        print("##ALT=<ID=BND,Description=\"Breakend\">", file=vcf_output)
    print("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">", file=vcf_output)
    print("##INFO=<ID=CUTPASTE,Number=0,Type=Flag,Description=\"Genomic origin of interspersed duplication seems to be deleted\">", file=vcf_output)
    print("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">", file=vcf_output)
    print("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">", file=vcf_output)    
    if options.read_names:
        print("##INFO=<ID=READS,Number=.,Type=String,Description=\"Names of all supporting reads\">", file=vcf_output)
    print("##FILTER=<ID=not_fully_covered,Description=\"Tandem duplication is not fully covered by a contig\">", file=vcf_output)
    print("##FILTER=<ID=incomplete_inversion,Description=\"Only one inversion breakpoint is supported\">", file=vcf_output)
    print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", file=vcf_output)
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + options.sample, file=vcf_output)


    # Prepare VCF entries depending on command-line parameters
    vcf_entries = []
    sequence_alleles = not options.symbolic_alleles
    if "DEL" in types_to_output:
        for candidate in deletion_candidates:
            vcf_entries.append((candidate.get_source(), candidate.get_vcf_entry(sequence_alleles, reference, options.read_names), "DEL"))
    if "INV" in types_to_output:
        for candidate in inversion_candidates:
            vcf_entries.append((candidate.get_source(), candidate.get_vcf_entry(sequence_alleles, reference, options.read_names), "INV"))
    if "INS" in types_to_output:
        for candidate in insertion_candidates:
            vcf_entries.append((candidate.get_destination(), candidate.get_vcf_entry(sequence_alleles, reference, options.read_names), "INS"))
    if options.duplications_as_insertions:
        if "INS" in types_to_output:
            for candidate in tandem_duplication_candidates:
                vcf_entries.append((candidate.get_destination(), candidate.get_vcf_entry_as_ins(sequence_alleles, reference, options.read_names), "INS"))
            for candidate in int_duplication_candidates:
                vcf_entries.append((candidate.get_destination(), candidate.get_vcf_entry_as_ins(sequence_alleles, reference, options.read_names), "INS"))
    else:
        if "DUP_TAN" in types_to_output:
            for candidate in tandem_duplication_candidates:
                vcf_entries.append((candidate.get_source(), candidate.get_vcf_entry_as_dup(options.read_names), "DUP_TAN"))
        if "DUP_INT" in types_to_output:
            for candidate in int_duplication_candidates:
                vcf_entries.append((candidate.get_source(), candidate.get_vcf_entry_as_dup(options.read_names), "DUP_INT"))
    if "BND" in types_to_output:
        for candidate in breakend_candidates:
            vcf_entries.append(((candidate.get_source()[0], candidate.get_source()[1], candidate.get_source()[1] + 1), candidate.get_vcf_entry(options.read_names), "BND"))

    if sequence_alleles:
        reference.close()

    # Sort and write entries to VCF
    svtype_counter = defaultdict(int)
    for source, entry, svtype in sorted_nicely(vcf_entries):
        variant_id = "svim_asm.{svtype}.{number}".format(svtype = svtype, number = svtype_counter[svtype] + 1)
        entry_with_id = entry.replace("PLACEHOLDERFORID", variant_id, 1)
        svtype_counter[svtype] += 1
        print(entry_with_id, file=vcf_output)

    vcf_output.close()