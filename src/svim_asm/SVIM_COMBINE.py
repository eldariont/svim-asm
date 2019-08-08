import os
import logging
import re

from collections import defaultdict
from math import pow, sqrt
import time
from statistics import mean, stdev
from pysam import FastaFile

from svim_asm.SVCandidate import CandidateInversion, CandidateDuplicationTandem, CandidateDeletion, CandidateInsertion, CandidateBreakend


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
                    options):
    vcf_output = open(options.working_dir + '/final_results.vcf', 'w')

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
    print("##FILTER=<ID=not_fully_covered,Description=\"Tandem duplication is not fully covered by a single read\">", file=vcf_output)
    print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", file=vcf_output)
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + options.sample, file=vcf_output)

    # Open reference genome sequence file
    sequence_alleles = not options.symbolic_alleles
    if sequence_alleles:
        try:
            reference = FastaFile(options.genome)
        except ValueError:
            logging.warning("The given reference genome is missing an index file ({path}.fai). Sequence alleles cannot be retrieved.".format(options.genome))
            sequence_alleles = False
        except IOError:
            logging.warning("The given reference genome is missing ({path}). Sequence alleles cannot be retrieved.".format(options.genome))
            sequence_alleles = False
    else:
        reference = None

    # Prepare VCF entries depending on command-line parameters
    vcf_entries = []
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