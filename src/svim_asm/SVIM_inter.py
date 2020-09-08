from __future__ import print_function

import sys
from statistics import mean
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster


from svim_asm.SVCandidate import CandidateDeletion, CandidateInsertion, CandidateInversion, CandidateBreakend, CandidateDuplicationTandem, CandidateDuplicationInterspersed


def is_similar(chr1, start1, end1, chr2, start2, end2):
    if chr1 == chr2 and abs(start1 - start2) < 20 and abs(end1 - end2) < 20:
        return True
    else:
        return False


def reciprocal_overlap_distance(inversion1, inversion2):
    start1, end1, direction1 = inversion1
    start2, end2, direction2 = inversion2
    #Inversion breakpoints with same direction cannot be joined
    if direction1 == direction2:
        return 1
    #Non-overlapping inversion breakpoints cannot be joined
    if start2 >= end1:
        return 1
    if start1 >= end2:
        return 1
    
    if start2 >= start1:
        overlap = min(end1, end2) - start2
    else:
        overlap = min(end1, end2) - start1
    
    relative_overlap1 = overlap / float(end1 - start1)
    relative_overlap2 = overlap / float(end2 - start2)
    minimum_relative_overlap = min(relative_overlap1, relative_overlap2)
    return 1 - minimum_relative_overlap


def process_overlapping_inversions(active_inversions, query_name, bam):
    if len(active_inversions) < 2:
        clusters = [active_inversions]
    else:
        data = np.array( [[inversion[1], inversion[2], 0 if inversion[3].split("_")[0] == "left" else 1] for inversion in active_inversions])
        Z = linkage(data, method = "complete", metric = reciprocal_overlap_distance)
        cluster_indices = list(fcluster(Z, 0.3, criterion='distance'))
        clusters = [[] for i in range(max(cluster_indices))]
        for inversion_index, cluster_index in enumerate(cluster_indices):
            clusters[cluster_index-1].append(active_inversions[inversion_index])

    inversion_candidates = []
    for cluster in clusters:
        chrom = cluster[0][0]
        start = max([i[1] for i in cluster])
        end = min([i[2] for i in cluster])
        complete = True if len(cluster) > 1 else False
        inversion_candidates.append(CandidateInversion(chrom, start, end, [query_name], complete, bam))
    return inversion_candidates

def analyze_read_segments(primary, supplementaries, bam, options):
    read_name = primary.query_name
    alignments = [primary] + supplementaries
    alignment_list = []
    for alignment in alignments:
        #correct query coordinates for reversely mapped reads
        if alignment.is_reverse:
            q_start = alignment.infer_read_length() - alignment.query_alignment_end
            q_end = alignment.infer_read_length() - alignment.query_alignment_start
        else:
            q_start = alignment.query_alignment_start
            q_end = alignment.query_alignment_end

        new_alignment_dict = {  'q_start': q_start, 
                                'q_end': q_end, 
                                'ref_id': alignment.reference_id, 
                                'ref_start': alignment.reference_start, 
                                'ref_end': alignment.reference_end,
                                'is_reverse': alignment.is_reverse  }
        alignment_list.append(new_alignment_dict)

    sorted_alignment_list = sorted(alignment_list, key=lambda aln: (aln['q_start'], aln['q_end']))
    #inferred_read_length = alignments[0].infer_read_length()

    sv_candidates = []
    tandem_duplications = []
    translocations = []
    inversions = []

    for index in range(len(sorted_alignment_list) - 1):
        alignment_current = sorted_alignment_list[index]
        alignment_next = sorted_alignment_list[index + 1]

        distance_on_read = alignment_next['q_start'] - alignment_current['q_end']

        #Same chromosome
        if alignment_current['ref_id'] == alignment_next['ref_id']:
            ref_chr = bam.get_reference_name(alignment_current['ref_id'])
            #Same orientation
            if alignment_current['is_reverse'] == alignment_next['is_reverse']:
                #Compute distance on reference depending on orientation
                if alignment_current['is_reverse']:
                    distance_on_reference = alignment_current['ref_start'] - alignment_next['ref_end']
                else:
                    distance_on_reference = alignment_next['ref_start'] - alignment_current['ref_end']
                #No overlap on read
                if distance_on_read >= -options.query_overlap_tolerance:
                    #No overlap on reference
                    if distance_on_reference >= -options.reference_overlap_tolerance:
                        deviation = distance_on_read - distance_on_reference
                        #INS candidate
                        if deviation >= options.min_sv_size:
                            #No gap on reference
                            if distance_on_reference <= options.reference_gap_tolerance:
                                if not alignment_current['is_reverse']:
                                    insertion_seq = primary.query_sequence[alignment_current['q_end']:alignment_current['q_end']+deviation]
                                    sv_candidates.append(CandidateInsertion(ref_chr, alignment_current['ref_end'], alignment_current['ref_end'] + deviation, [read_name], insertion_seq, bam))
                                else:
                                    insertion_seq = primary.query_sequence[primary.infer_read_length() - alignment_next['q_start']:primary.infer_read_length() - alignment_next['q_start'] + deviation]
                                    sv_candidates.append(CandidateInsertion(ref_chr, alignment_current['ref_start'], alignment_current['ref_start'] + deviation, [read_name], insertion_seq, bam))
                        #DEL candidate
                        elif -options.max_sv_size <= deviation <= -options.min_sv_size:
                            #No gap on read
                            if distance_on_read <= options.query_gap_tolerance:
                                if not alignment_current['is_reverse']:
                                    sv_candidates.append(CandidateDeletion(ref_chr, alignment_current['ref_end'], alignment_current['ref_end'] - deviation, [read_name], bam))
                                else:
                                    sv_candidates.append(CandidateDeletion(ref_chr, alignment_next['ref_end'], alignment_next['ref_end'] - deviation, [read_name], bam))
                        #Either very large DEL or TRANS
                        elif deviation < -options.max_sv_size:
                            #No gap on read
                            if distance_on_read <= options.query_gap_tolerance:
                                if not alignment_current['is_reverse']:
                                    sv_candidates.append(CandidateBreakend(ref_chr, alignment_current['ref_end'] - 1, 'fwd', ref_chr, alignment_next['ref_start'], 'fwd', [read_name], bam))
                                    translocations.append(('fwd', 'fwd', ref_chr, alignment_current['ref_end'] - 1, ref_chr, alignment_next['ref_start']))
                                else:
                                    sv_candidates.append(CandidateBreakend(ref_chr, alignment_current['ref_start'], 'rev', ref_chr, alignment_next['ref_end'] - 1, 'rev', [read_name], bam))
                                    translocations.append(('rev', 'rev', ref_chr, alignment_current['ref_start'], ref_chr, alignment_next['ref_end'] - 1))
                    #overlap on reference
                    else:
                        #No gap on read
                        if distance_on_read <= options.query_gap_tolerance:
                            deviation = distance_on_read - distance_on_reference
                            #Tandem Duplication
                            if deviation >= options.min_sv_size:
                                if not alignment_current['is_reverse']:
                                    #Tandem Duplication (fully covered)
                                    if alignment_next['ref_end'] > alignment_current['ref_start']:
                                        tandem_duplications.append((ref_chr, alignment_next['ref_start'], alignment_next['ref_start'] + deviation, True, True))
                                    #Tandem duplication (not fully covered)
                                    elif distance_on_reference >= -options.max_sv_size:
                                        tandem_duplications.append((ref_chr, alignment_next['ref_start'], alignment_next['ref_start'] + deviation, False, True))
                                    #Either very large TANDEM or TRANS
                                    else:
                                        sv_candidates.append(CandidateBreakend(ref_chr, alignment_current['ref_end'] - 1, 'fwd', ref_chr, alignment_next['ref_start'], 'fwd', [read_name], bam))
                                        translocations.append(('fwd', 'fwd', ref_chr, alignment_current['ref_end'] - 1, ref_chr, alignment_next['ref_start']))
                                else:
                                    #Tandem Duplication
                                    if alignment_next['ref_start'] < alignment_current['ref_end']:
                                        tandem_duplications.append((ref_chr, alignment_current['ref_start'], alignment_current['ref_start'] + deviation, True, False))
                                    #Large tandem duplication
                                    elif distance_on_reference >= -options.max_sv_size:
                                        tandem_duplications.append((ref_chr, alignment_current['ref_start'], alignment_current['ref_start'] + deviation, False, False))
                                    #Either very large TANDEM or TRANS
                                    else:
                                        sv_candidates.append(CandidateBreakend(ref_chr, alignment_current['ref_start'], 'rev', ref_chr, alignment_next['ref_end'] - 1, 'rev', [read_name], bam))
                                        translocations.append(('rev', 'rev', ref_chr, alignment_current['ref_start'], ref_chr, alignment_next['ref_end'] - 1))
            #Different orientations
            else:
                #Normal to reverse
                if not alignment_current['is_reverse'] and alignment_next['is_reverse']:
                    distance_on_reference = alignment_next['ref_end'] - alignment_current['ref_end']
                    deviation = distance_on_read - distance_on_reference
                    if -options.query_overlap_tolerance <= distance_on_read <= options.query_gap_tolerance:
                        if alignment_next['ref_start'] - alignment_current['ref_end'] >= -options.reference_overlap_tolerance: # Case 1
                            #INV candidate
                            if options.min_sv_size <= -deviation <= options.max_sv_size:
                                inversions.append((ref_chr, alignment_current['ref_end'], alignment_current['ref_end'] - deviation, "left_fwd"))
                                #transitions.append(('inversion', 'left_fwd', ref_chr, alignment_current['ref_end'], alignment_next['ref_end']))
                            #Either very large INV or TRANS
                            else:
                                sv_candidates.append(CandidateBreakend(ref_chr, alignment_current['ref_end'] - 1, 'fwd', ref_chr, alignment_next['ref_end'] - 1, 'rev', [read_name], bam))
                                translocations.append(('fwd', 'rev', ref_chr, alignment_current['ref_end'] - 1, ref_chr, alignment_next['ref_end'] - 1))
                        elif alignment_current['ref_start'] - alignment_next['ref_end'] >= -options.reference_overlap_tolerance: # Case 3
                            #INV candidate
                            if options.min_sv_size <= deviation <= options.max_sv_size:
                                inversions.append((ref_chr, alignment_next['ref_end'], alignment_next['ref_end'] + deviation, "left_rev"))
                                #transitions.append(('inversion', 'left_rev', ref_chr, alignment_next['ref_end'], alignment_current['ref_end']))
                            #Either very large INV or TRANS
                            else:
                                sv_candidates.append(CandidateBreakend(ref_chr, alignment_current['ref_end'] - 1, 'fwd', ref_chr, alignment_next['ref_end'] - 1, 'rev', [read_name], bam))
                                translocations.append(('fwd', 'rev', ref_chr, alignment_current['ref_end'] - 1, ref_chr, alignment_next['ref_end'] - 1))
                    else:
                        pass
                        #print("Overlapping read segments in read", read_name)
                #Reverse to normal
                if alignment_current['is_reverse'] and not alignment_next['is_reverse']:
                    distance_on_reference = alignment_next['ref_start'] - alignment_current['ref_start'] 
                    deviation = distance_on_read - distance_on_reference
                    if -options.query_overlap_tolerance <= distance_on_read <= options.query_gap_tolerance:
                        if alignment_next['ref_start'] - alignment_current['ref_end'] >= -options.reference_overlap_tolerance: # Case 2
                            #INV candidate
                            if options.min_sv_size <= -deviation <= options.max_sv_size:
                                inversions.append((ref_chr, alignment_current['ref_start'], alignment_current['ref_start'] - deviation, "right_fwd"))
                                #transitions.append(('inversion', 'right_fwd', ref_chr, alignment_current['ref_start'], alignment_next['ref_start']))
                            #Either very large INV or TRANS
                            else:
                                sv_candidates.append(CandidateBreakend(ref_chr, alignment_current['ref_start'], 'rev', ref_chr, alignment_next['ref_start'], 'fwd', [read_name], bam))
                                translocations.append(('rev', 'fwd', ref_chr, alignment_current['ref_start'], ref_chr, alignment_next['ref_start']))
                        elif alignment_current['ref_start'] - alignment_next['ref_end'] >= -options.reference_overlap_tolerance: # Case 4
                            #INV candidate
                            if options.min_sv_size <= deviation <= options.max_sv_size:
                                inversions.append((ref_chr, alignment_next['ref_start'], alignment_next['ref_start'] + deviation, "right_rev"))
                                #transitions.append(('inversion', 'right_rev', ref_chr, alignment_next['ref_start'], alignment_current['ref_start']))
                            #Either very large INV or TRANS
                            else:
                                sv_candidates.append(CandidateBreakend(ref_chr, alignment_current['ref_start'], 'rev', ref_chr, alignment_next['ref_start'], 'fwd', [read_name], bam))
                                translocations.append(('rev', 'fwd', ref_chr, alignment_current['ref_start'], ref_chr, alignment_next['ref_start']))
                    else:
                        pass
                        #print("Overlapping read segments in read", read_name)
        #Different chromosomes
        else:
            ref_chr_current = bam.getrname(alignment_current['ref_id'])
            ref_chr_next = bam.getrname(alignment_next['ref_id'])
            #Same orientation
            if alignment_current['is_reverse'] == alignment_next['is_reverse']:
                #No overlap on read
                if distance_on_read >= -options.query_overlap_tolerance:
                    #No gap on read
                    if distance_on_read <= options.query_gap_tolerance:
                        if not alignment_current['is_reverse']:
                            sv_candidates.append(CandidateBreakend(ref_chr_current, alignment_current['ref_end'] - 1, 'fwd', ref_chr_next, alignment_next['ref_start'], 'fwd', [read_name], bam))
                            translocations.append(('fwd', 'fwd', ref_chr_current, alignment_current['ref_end'] - 1, ref_chr_next, alignment_next['ref_start']))
                        else:
                            sv_candidates.append(CandidateBreakend(ref_chr_current, alignment_current['ref_start'], 'rev', ref_chr_next, alignment_next['ref_end'] - 1, 'rev', [read_name], bam))
                            translocations.append(('rev', 'rev', ref_chr_current, alignment_current['ref_start'], ref_chr_next, alignment_next['ref_end'] - 1))
                #Overlap on read
                else:
                    pass
                    #print("Overlapping read segments in read", read_name)
            #Different orientation
            else:
                #No overlap on read
                if distance_on_read >= -options.query_overlap_tolerance:
                    #No gap on read
                    if distance_on_read <= options.query_gap_tolerance:
                        if not alignment_current['is_reverse']:
                            sv_candidates.append(CandidateBreakend(ref_chr_current, alignment_current['ref_end'] - 1, 'fwd', ref_chr_next, alignment_next['ref_end'] - 1, 'rev', [read_name], bam))
                            translocations.append(('fwd', 'rev', ref_chr_current, alignment_current['ref_end'] - 1, ref_chr_next, alignment_next['ref_end'] - 1))
                        else:
                            sv_candidates.append(CandidateBreakend(ref_chr_current, alignment_current['ref_start'], 'rev', ref_chr_next, alignment_next['ref_start'], 'fwd', [read_name], bam))
                            translocations.append(('rev', 'fwd', ref_chr_current, alignment_current['ref_start'], ref_chr_next, alignment_next['ref_start']))
                #Overlap on read
                else:
                    pass
                    #print("Overlapping read segments in read", read_name)

    #Handle tandem duplications
    current_chromosome = None
    current_starts = []
    current_ends = []
    current_copy_number = 0
    current_fully_covered = []
    for tandem_duplication in tandem_duplications:
        if current_chromosome == None:
            current_chromosome = tandem_duplication[0]
            current_starts.append(tandem_duplication[1])
            current_ends.append(tandem_duplication[2])
            current_copy_number = 1
            current_fully_covered.append(tandem_duplication[3])
            current_direction = tandem_duplication[4]
        else:
            if is_similar(current_chromosome, mean(current_starts), mean(current_ends), tandem_duplication[0], tandem_duplication[1], tandem_duplication[2]) and current_direction == tandem_duplication[4]:
                current_starts.append(tandem_duplication[1])
                current_ends.append(tandem_duplication[2])
                current_copy_number += 1
                current_fully_covered.append(tandem_duplication[3])
            else:
                fully_covered = True if sum(current_fully_covered) else False
                sv_candidates.append(CandidateDuplicationTandem(current_chromosome, int(mean(current_starts)), int(mean(current_ends)), current_copy_number, fully_covered, [read_name], bam))
                current_chromosome = tandem_duplication[0]
                current_starts =[tandem_duplication[1]]
                current_ends =[tandem_duplication[2]]
                current_copy_number = 1
                current_fully_covered = [tandem_duplication[3]]
    if current_chromosome != None:
        fully_covered = True if sum(current_fully_covered) else False
        sv_candidates.append(CandidateDuplicationTandem(current_chromosome, int(mean(current_starts)), int(mean(current_ends)), current_copy_number, fully_covered, [read_name], bam))

    #Handle interspersed duplications
    for this_index in range(len(translocations)):
        this_dir1 = translocations[this_index][0]
        this_dir2 = translocations[this_index][1]
        this_chr1 = translocations[this_index][2]
        this_pos1 = translocations[this_index][3]
        this_chr2 = translocations[this_index][4]
        this_pos2 = translocations[this_index][5]

        for before_dir1, before_dir2, before_chr1, before_pos1, before_chr2, before_pos2 in translocations[:this_index]:
            #Same direction at destination and origin
            if before_dir1 == this_dir2 and before_dir2 == this_dir1:
                #Same position at destination
                if is_similar(before_chr1, before_pos1, 0, this_chr2, this_pos2, 0):
                    #Same chromosome for origin
                    if before_chr2 == this_chr1:
                        #INS_DUP candidate
                        if before_dir2 == before_dir1:
                            if before_dir1 == 'fwd':
                                length = this_pos1 + 1 - before_pos2
                                if options.min_sv_size <= length <= options.max_sv_size:
                                    sv_candidates.append(CandidateDuplicationInterspersed(before_chr2, before_pos2, this_pos1 + 1, before_chr1, int(mean([before_pos1 + 1, this_pos2])), int(mean([before_pos1 + 1, this_pos2])) + length, [read_name], bam))
                            elif before_dir1 == 'rev':
                                length = before_pos2 + 1 - this_pos1
                                if options.min_sv_size <= length <= options.max_sv_size:
                                    sv_candidates.append(CandidateDuplicationInterspersed(before_chr2, this_pos1, before_pos2 + 1, before_chr1, int(mean([before_pos1, this_pos2 + 1])), int(mean([before_pos1, this_pos2 + 1])) + length, [read_name], bam))
                        #INV_INS_DUP candidate
                        else:
                            pass

    #Handle inversions (simple inversions produce two novel adjacencies that need to be merged for a complete candidate)
    sorted_inversions = sorted(inversions, key=lambda inversion: (inversion[0], inversion[1], inversion[2]))
    active_inversions = []
    for inversion in sorted_inversions:
        chrom, start, end, direction = inversion
        if len(active_inversions) == 0:
            active_inversions.append(inversion)
        else:
            #If current inversion overlaps one of the active inversions
            if chrom == active_inversions[-1][0] and start < max([i[2] for i in active_inversions]):
                active_inversions.append(inversion)
            else:
                #Cluster inversions
                sv_candidates.extend(process_overlapping_inversions(active_inversions, read_name, bam))
                active_inversions = []
    if len(active_inversions) > 0:
        sv_candidates.extend(process_overlapping_inversions(active_inversions, read_name, bam))  

    return sv_candidates
