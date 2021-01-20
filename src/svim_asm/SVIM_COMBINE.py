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
        chr_length = reference.get_reference_length(region_chr)
        region_start = max(0, min(candidate1.source_start, candidate2.source_start) - 100)
        region_end = min(chr_length, max(candidate1.source_end, candidate2.source_end) + 100)
        haplotype1 = reference.fetch(region_chr, region_start, candidate1.source_start).upper() + reference.fetch(region_chr, candidate1.source_end, region_end).upper()
        haplotype2 = reference.fetch(region_chr, region_start, candidate2.source_start).upper() + reference.fetch(region_chr, candidate2.source_end, region_end).upper()
        editDistance = align(haplotype1, haplotype2)["editDistance"]
    elif candidate1.type == "INV":
        region_chr = candidate1.source_contig
        chr_length = reference.get_reference_length(region_chr)
        region_start = max(0, min(candidate1.source_start, candidate2.source_start) - 100)
        region_end = min(chr_length, max(candidate1.source_end, candidate2.source_end) + 100)
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
        chr_length = reference.get_reference_length(region_chr)
        region_start = max(0, min(candidate1.dest_start, candidate2.dest_start) - 100)
        region_end = min(chr_length, max(candidate1.dest_start, candidate2.dest_start) + 100)
        haplotype1 = reference.fetch(region_chr, region_start, candidate1.dest_start).upper() + \
                     candidate1.sequence + \
                     reference.fetch(region_chr, candidate1.dest_start, region_end).upper()
        haplotype2 = reference.fetch(region_chr, region_start, candidate2.dest_start).upper() + \
                     candidate2.sequence + \
                     reference.fetch(region_chr, candidate2.dest_start, region_end).upper()
        editDistance = align(haplotype1, haplotype2)["editDistance"]
    elif candidate1.type == "DUP_TAN":
        region_chr = candidate1.source_contig
        chr_length = reference.get_reference_length(region_chr)
        region_start = max(0, min(candidate1.source_start, candidate2.source_start) - 100)
        region_end = min(chr_length, max(candidate1.source_end, candidate2.source_end) + 100)
        haplotype1 = reference.fetch(region_chr, region_start, candidate1.source_start).upper() + \
                     reference.fetch(region_chr, candidate1.source_start, candidate1.source_end).upper() * (candidate1.copies + 1) + \
                     reference.fetch(region_chr, candidate1.source_end, region_end).upper()
        haplotype2 = reference.fetch(region_chr, region_start, candidate2.source_start).upper() + \
                     reference.fetch(region_chr, candidate2.source_start, candidate2.source_end).upper() * (candidate2.copies + 1) + \
                     reference.fetch(region_chr, candidate2.source_end, region_end).upper()
        editDistance = align(haplotype1, haplotype2)["editDistance"]
    elif candidate1.type == "DUP_INT":
        region_chr = candidate1.dest_contig
        chr_length = reference.get_reference_length(region_chr)
        region_start = max(0, min(candidate1.dest_start, candidate2.dest_start) - 100)
        region_end = min(chr_length, max(candidate1.dest_start, candidate2.dest_start) + 100)
        haplotype1 = reference.fetch(region_chr, region_start, candidate1.dest_start).upper() + \
                     reference.fetch(candidate1.source_contig, candidate1.source_start, candidate1.source_end).upper() + \
                     reference.fetch(region_chr, candidate1.dest_start, region_end).upper()
        haplotype2 = reference.fetch(region_chr, region_start, candidate2.dest_start).upper() + \
                     reference.fetch(candidate2.source_contig, candidate2.source_start, candidate2.source_end).upper() + \
                     reference.fetch(region_chr, candidate2.dest_start, region_end).upper()
        editDistance = align(haplotype1, haplotype2)["editDistance"]

    return editDistance


def span_position_distance_breakends(candidate1, candidate2):
    candidate1_hap, candidate1_pos1, candidate1_dir1, candidate1_pos2, candidate1_dir2 = candidate1
    candidate2_hap, candidate2_pos1, candidate2_dir1, candidate2_pos2, candidate2_dir2 = candidate2
    if candidate1_hap != candidate2_hap:
        if candidate1_dir1 == candidate2_dir1 and candidate1_dir2 == candidate2_dir2:
            dist1 = abs(candidate1_pos1 - candidate2_pos1)
            dist2 = abs(candidate1_pos2 - candidate2_pos2)
            position_distance = (dist1 + dist2) / 3000
        else:
            position_distance = 99999
    else:
        position_distance = 99999
    return position_distance


def pair_haplotypes(partitions, reference, edit_distance_threshold = 10):
    clusters_final = []
    for partition in partitions:
        if len(partition) < 2:
            new_clusters = [partition]
        #Ignore very large partitions because they tend to be in difficult regions
        elif len(partition) > 10:
            logging.debug("Ignored partition of size {0} and type {1}: {2}".format(len(partition), partition[0][1].get_key()[0], ",".join(["{0}:{1}".format(sig[1].get_key()[1], sig[1].get_key()[2]) for sig in partition])))
            continue
        else:
            distances = []
            for i in range(len(partition)-1):
                for j in range(i+1, len(partition)):
                    distances.append(compute_distance(partition[i], partition[j], reference))
            Z = linkage(np.array(distances), method = "complete")
            cluster_indices = list(fcluster(Z, edit_distance_threshold, criterion='distance'))
            new_clusters = [[] for i in range(max(cluster_indices))]
            for candidate_index, cluster_index in enumerate(cluster_indices):
                new_clusters[cluster_index-1].append(partition[candidate_index])
        clusters_final.extend(new_clusters)
    return clusters_final


def pair_haplotypes_breakends(partitions, span_position_distance_threshold = 0.3):
    """Finds clusters in partitions using span-position distance and hierarchical clustering. 
    Assumes that all signatures in the given partition are of the same type and on the same contig"""
    clusters_final = []
    for partition in partitions:
        if len(partition) < 2:
            new_clusters = [partition]
        #Ignore very large partitions because they tend to be in difficult regions
        elif len(partition) > 10:
            continue
        else:
            data = np.array( [[haplotype, candidate.get_source()[1], 1 if candidate.source_direction == 'fwd' else 0, candidate.get_destination()[1], 1 if candidate.dest_direction == 'fwd' else 0] for (haplotype, candidate) in partition])
            Z = linkage(data, method = "complete", metric = span_position_distance_breakends)
            cluster_indices = list(fcluster(Z, span_position_distance_threshold, criterion='distance'))
            new_clusters = [[] for i in range(max(cluster_indices))]
            for candidate_index, cluster_index in enumerate(cluster_indices):
                new_clusters[cluster_index-1].append(partition[candidate_index])
        clusters_final.extend(new_clusters)
    return clusters_final


def pair_candidates(sv_candidates1, sv_candidates2, reference, bam, options):
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
    partitions = form_partitions(deletion_candidates1 + deletion_candidates2, options.partition_max_distance)
    clusters = pair_haplotypes(partitions, reference, options.max_edit_distance)
    for cluster in clusters:
        if len(cluster) == 1:
            candidate = cluster[0][1]
            genotype = "1/0" if cluster[0][0] == 1 else "0/1"
            paired_candidates.append(CandidateDeletion(candidate.source_contig, 
                                                       candidate.source_start, 
                                                       candidate.source_end, 
                                                       candidate.reads,
                                                       bam,
                                                       genotype))
        elif len(cluster) == 2:
            candidate = cluster[0][1]
            reads = cluster[0][1].reads + cluster[1][1].reads
            genotype = "1/1"
            paired_candidates.append(CandidateDeletion(candidate.source_contig, 
                                                       candidate.source_start, 
                                                       candidate.source_end, 
                                                       reads,
                                                       bam,
                                                       genotype))
        else:
            logging.error("Cluster size should be either 1 or 2 but is " + str(len(cluster)))
    
    #INVERSIONS
    logging.info("Pairing {0} inversions...".format(len(inversion_candidates1) + len(inversion_candidates2)))
    partitions = form_partitions(inversion_candidates1 + inversion_candidates2, options.partition_max_distance)
    clusters = pair_haplotypes(partitions, reference, options.max_edit_distance)
    for cluster in clusters:
        if len(cluster) == 1:
            candidate = cluster[0][1]
            genotype = "1/0" if cluster[0][0] == 1 else "0/1"
            paired_candidates.append(CandidateInversion(candidate.source_contig, 
                                                        candidate.source_start, 
                                                        candidate.source_end, 
                                                        candidate.reads, 
                                                        candidate.complete, 
                                                        bam, 
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
                                                        bam, 
                                                        genotype))
        else:
            logging.error("Cluster size should be either 1 or 2 but is " + str(len(cluster)))

    #INSERTIONS
    logging.info("Pairing {0} insertions...".format(len(insertion_candidates1) + len(insertion_candidates2)))
    partitions = form_partitions(insertion_candidates1 + insertion_candidates2, options.partition_max_distance)
    clusters = pair_haplotypes(partitions, reference, options.max_edit_distance)
    for cluster in clusters:
        if len(cluster) == 1:
            candidate = cluster[0][1]
            genotype = "1/0" if cluster[0][0] == 1 else "0/1"
            paired_candidates.append(CandidateInsertion(candidate.dest_contig, 
                                                        candidate.dest_start, 
                                                        candidate.dest_end, 
                                                        candidate.reads, 
                                                        candidate.sequence, 
                                                        bam, 
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
                                                        bam, 
                                                        genotype))
        else:
            logging.error("Cluster size should be either 1 or 2 but is " + str(len(cluster)))

    #TANDEM DUPLICATIONS
    logging.info("Pairing {0} tandem duplications...".format(len(tandem_duplication_candidates1) + len(tandem_duplication_candidates2)))
    partitions = form_partitions(tandem_duplication_candidates1 + tandem_duplication_candidates2, options.partition_max_distance)
    clusters = pair_haplotypes(partitions, reference, options.max_edit_distance)
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
                                                                bam, 
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
                                                                bam, 
                                                                genotype))
        else:
            logging.error("Cluster size should be either 1 or 2 but is " + str(len(cluster)))

    #INTERSPERSED DUPLICATIONS
    logging.info("Pairing {0} interspersed duplications...".format(len(interspersed_duplication_candidates1) + len(interspersed_duplication_candidates2)))
    partitions = form_partitions(interspersed_duplication_candidates1 + interspersed_duplication_candidates2, options.partition_max_distance)
    clusters = pair_haplotypes(partitions, reference, options.max_edit_distance)
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
                                                                      bam, 
                                                                      candidate.cutpaste,
                                                                      genotype))
        elif len(cluster) == 2:
            candidate = cluster[0][1]
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
                                                                      bam, 
                                                                      cutpaste,
                                                                      genotype))
        else:
            logging.error("Cluster size should be either 1 or 2 but is " + str(len(cluster)))

    #BREAKENDS
    logging.info("Pairing {0} breakends...".format(len(breakend_candidates1) + len(breakend_candidates2)))
    partitions = form_partitions(breakend_candidates1 + breakend_candidates2, options.partition_max_distance)
    clusters = pair_haplotypes_breakends(partitions)
    for cluster in clusters:
        if len(cluster) == 1:
            candidate = cluster[0][1]
            genotype = "1/0" if cluster[0][0] == 1 else "0/1"
            paired_candidates.append(CandidateBreakend(candidate.source_contig, 
                                                        candidate.source_start, 
                                                        candidate.source_direction, 
                                                        candidate.dest_contig, 
                                                        candidate.dest_start, 
                                                        candidate.dest_direction,
                                                        candidate.reads,
                                                        bam,
                                                        genotype))
        elif len(cluster) == 2:
            candidate = cluster[0][1]
            reads = cluster[0][1].reads + cluster[1][1].reads
            genotype = "1/1"
            paired_candidates.append(CandidateBreakend(candidate.source_contig, 
                                                        candidate.source_start, 
                                                        candidate.source_direction, 
                                                        candidate.dest_contig, 
                                                        candidate.dest_start, 
                                                        candidate.dest_direction,
                                                        reads,
                                                        bam,
                                                        genotype))
        else:
            logging.error("Cluster size should be either 1 or 2 but is " + str(len(cluster)))
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
    if (not options.tandem_duplications_as_insertions and "DUP:TANDEM" in types_to_output) or \
       (not options.interspersed_duplications_as_insertions and "DUP:INT" in types_to_output):
        print("##ALT=<ID=DUP,Description=\"Duplication\">", file=vcf_output)
    if not options.tandem_duplications_as_insertions and "DUP:TANDEM" in types_to_output:
        print("##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">", file=vcf_output)
    if not options.interspersed_duplications_as_insertions and "DUP:INT" in types_to_output:
        print("##ALT=<ID=DUP:INT,Description=\"Interspersed Duplication\">", file=vcf_output)
    if "INS" in types_to_output:
        print("##ALT=<ID=INS,Description=\"Insertion\">", file=vcf_output)
    if "BND" in types_to_output:
        print("##ALT=<ID=BND,Description=\"Breakend\">", file=vcf_output)
    print("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">", file=vcf_output)
    print("##INFO=<ID=CUTPASTE,Number=0,Type=Flag,Description=\"Genomic origin of interspersed duplication seems to be deleted\">", file=vcf_output)
    print("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">", file=vcf_output)
    print("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">", file=vcf_output)    
    if options.query_names:
        print("##INFO=<ID=READS,Number=.,Type=String,Description=\"Names of all supporting reads\">", file=vcf_output)
    print("##FILTER=<ID=not_fully_covered,Description=\"Tandem duplication is not fully covered by a contig\">", file=vcf_output)
    print("##FILTER=<ID=incomplete_inversion,Description=\"Only one inversion breakpoint is supported\">", file=vcf_output)
    print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", file=vcf_output)
    if not options.tandem_duplications_as_insertions and "DUP:TANDEM" in types_to_output:
        print("##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number of tandem duplication (e.g. 2 for one additional copy)\">", file=vcf_output)
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + options.sample, file=vcf_output)


    # Prepare VCF entries depending on command-line parameters
    vcf_entries = []
    sequence_alleles = not options.symbolic_alleles
    if "DEL" in types_to_output:
        for candidate in deletion_candidates:
            contig, start, end = candidate.get_source()
            vcf_entries.append(((contig, max(1, start), end), candidate.get_vcf_entry(sequence_alleles, reference, options.query_names), "DEL"))
    if "INV" in types_to_output:
        for candidate in inversion_candidates:
            contig, start, end = candidate.get_source()
            vcf_entries.append(((contig, start+1, end), candidate.get_vcf_entry(sequence_alleles, reference, options.query_names), "INV"))
    if "INS" in types_to_output:
        for candidate in insertion_candidates:
            contig, start, end = candidate.get_destination()
            vcf_entries.append(((contig, max(1, start), end), candidate.get_vcf_entry(sequence_alleles, reference, options.query_names), "INS"))
    if options.tandem_duplications_as_insertions:
        if "INS" in types_to_output:
            for candidate in tandem_duplication_candidates:
                vcf_entries.append(((candidate.source_contig, candidate.source_start+1, candidate.source_end), candidate.get_vcf_entry_as_ins(sequence_alleles, reference, options.query_names), "INS"))
    else:
        if "DUP:TANDEM" in types_to_output:
            for candidate in tandem_duplication_candidates:
                vcf_entries.append(((candidate.source_contig, candidate.source_start+1, candidate.source_end), candidate.get_vcf_entry_as_dup(options.query_names), "DUP_TANDEM"))
    if options.interspersed_duplications_as_insertions:
        if "INS" in types_to_output:
            for candidate in int_duplication_candidates:
                contig, start, end = candidate.get_destination()
                vcf_entries.append(((contig, max(1, start), end), candidate.get_vcf_entry_as_ins(sequence_alleles, reference, options.query_names), "INS"))
    else:
        if "DUP:INT" in types_to_output:
            for candidate in int_duplication_candidates:
                contig, start, end = candidate.get_source()
                vcf_entries.append(((contig, start+1, end), candidate.get_vcf_entry_as_dup(options.query_names), "DUP_INT"))
    if "BND" in types_to_output:
        for candidate in breakend_candidates:
            vcf_entries.append(((candidate.get_source()[0], candidate.get_source()[1]+1, candidate.get_source()[1] + 2), candidate.get_vcf_entry(options.query_names), "BND"))
            vcf_entries.append(((candidate.get_destination()[0], candidate.get_destination()[1]+1, candidate.get_destination()[1] + 2), candidate.get_vcf_entry_reverse(options.query_names), "BND"))

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
