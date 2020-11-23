import logging
import pysam

from svim_asm.SVIM_intra import analyze_alignment_indel
from svim_asm.SVIM_inter import analyze_read_segments


def retrieve_other_alignments(main_alignment, bam):
    """Reconstruct other alignments of the same read for a given alignment from the SA tag"""
    #reconstructing other alignments from SA tag does not work if sequence of main_alignment is hard-clipped
    if main_alignment.get_cigar_stats()[0][5] > 0:
        return []
    try:
        sa_tag = main_alignment.get_tag("SA").split(";")
    except KeyError:
        return []
    other_alignments = []
    # For each other alignment encoded in the SA tag
    for element in sa_tag:
        # Read information from the tag
        fields = element.split(",")
        if len(fields) != 6:
            continue
        rname = fields[0]
        pos = int(fields[1])
        strand = fields[2]
        # CIGAR string encoded in SA tag is shortened
        cigar = fields[3]
        mapq = int(fields[4])
        nm = int(fields[5])

        # Generate an aligned segment from the information
        a = pysam.AlignedSegment()
        a.query_name = main_alignment.query_name
        a.query_sequence = ''
        if strand == "+":
            a.flag = 2048
        else:
            a.flag = 2064
        a.reference_id = bam.get_tid(rname)
        a.reference_start = pos - 1
        try:
            a.mapping_quality = mapq
        except OverflowError:
            a.mapping_quality = 0
        a.cigarstring = cigar
        a.next_reference_id = -1
        a.next_reference_start = -1
        a.template_length = 0
        a.query_qualities = main_alignment.query_qualities
        a.set_tags([("NM", nm, "i")])

        other_alignments.append(a)
    return other_alignments


def analyze_alignment_file_coordsorted(bam, options):
    chromosomes = bam.references
    sv_candidates = []
    for current_chromosome in chromosomes:
        alignment_it = bam.fetch(contig = current_chromosome)
        logging.info("Processing chromosome {0}...".format(current_chromosome))

        while True:
            try:
                current_alignment = next(alignment_it)
                if current_alignment.is_unmapped or current_alignment.is_secondary or current_alignment.mapping_quality < options.min_mapq:
                    continue
                if current_alignment.is_supplementary:
                    sv_candidates.extend(analyze_alignment_indel(current_alignment, bam, current_alignment.query_name, options))
                else:
                    supplementary_alignments = retrieve_other_alignments(current_alignment, bam)
                    good_suppl_alns = [aln for aln in supplementary_alignments if not aln.is_unmapped and aln.mapping_quality >= options.min_mapq]

                    sv_candidates.extend(analyze_alignment_indel(current_alignment, bam, current_alignment.query_name, options))
                    sv_candidates.extend(analyze_read_segments(current_alignment, good_suppl_alns, bam, options))
            except StopIteration:
                break
    return sv_candidates
