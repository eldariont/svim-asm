from __future__ import print_function

import sys
import os

from callPacCluster import partition_and_cluster_unilocal, partition_and_cluster_bilocal
from SVEvidence import EvidenceTranslocation
from callPacMerge import merge_insertions_from, merge_translocations_at_deletions, merge_translocations_at_insertions


def complete_translocations(translocation_evidences):
    """Generate a complete list of translocation by adding all reversed translocations"""

    reversed_translocations = []
    for evidence in translocation_evidences:
        reversed_translocations.append(EvidenceTranslocation(evidence.contig2, evidence.pos2, evidence.contig1, evidence.pos1, evidence.evidence, evidence.read))
    return translocation_evidences + reversed_translocations


def cluster_sv_evidences(sv_evidences):
    """Takes a list of SVEvidences and splits them up by type. The SVEvidences of each type are clustered and returned as a tuple of
    (deletion_evidence_clusters, insertion_evidence_clusters, inversion_evidence_clusters, tandem_duplication_evidence_clusters, insertion_from_evidence_clusters, completed_translocation_evidences)."""

    deletion_evidences = [ev for ev in sv_evidences if ev.type == 'del']
    insertion_evidences = [ev for ev in sv_evidences if ev.type == 'ins']
    inversion_evidences = [ev for ev in sv_evidences if ev.type == 'inv']
    tandem_duplication_evidences = [ev for ev in sv_evidences if ev.type == 'dup']
    translocation_evidences = [ev for ev in sv_evidences if ev.type == 'tra']
    insertion_from_evidences = [ev for ev in sv_evidences if ev.type == 'ins_dup']

    print("INFO: Found {0}/{1}/{2}/{3}/{4}/{5} evidences for deletions, insertions, inversions, tandem duplications, translocations, and insertion_from, respectively.".format(
        len(deletion_evidences), len(insertion_evidences), len(inversion_evidences), len(tandem_duplication_evidences), len(translocation_evidences), len(insertion_from_evidences)), file=sys.stderr)
    
    # Cluster SV evidences
    print("INFO: Cluster deletion evidences..", file=sys.stderr)
    deletion_evidence_clusters = partition_and_cluster_unilocal(deletion_evidences)
    print("INFO: Cluster insertion evidences..", file=sys.stderr)
    insertion_evidence_clusters = partition_and_cluster_unilocal(insertion_evidences)
    print("INFO: Cluster inversion evidences..", file=sys.stderr)
    inversion_evidence_clusters = partition_and_cluster_unilocal(inversion_evidences)
    print("INFO: Cluster tandem duplication evidences..", file=sys.stderr)
    tandem_duplication_evidence_clusters = partition_and_cluster_bilocal(tandem_duplication_evidences)
    print("INFO: Cluster insertion evidences with source..", file=sys.stderr)
    insertion_from_evidence_clusters = partition_and_cluster_bilocal(insertion_from_evidences)

    return (deletion_evidence_clusters, insertion_evidence_clusters, inversion_evidence_clusters, tandem_duplication_evidence_clusters, insertion_from_evidence_clusters, complete_translocations(translocation_evidences))


def write_evidence_clusters_bed(working_dir, clusters):
    """Write evidence clusters into working directory in BED format."""
    deletion_evidence_clusters, insertion_evidence_clusters, inversion_evidence_clusters, tandem_duplication_evidence_clusters, insertion_from_evidence_clusters, completed_translocations = clusters

    # Print SV evidence clusters
    if not os.path.exists(working_dir + '/evidences'):
        os.mkdir(working_dir + '/evidences')
    deletion_evidence_output = open(working_dir + '/evidences/del.bed', 'w')
    insertion_evidence_output = open(working_dir + '/evidences/ins.bed', 'w')
    inversion_evidence_output = open(working_dir + '/evidences/inv.bed', 'w')
    tandem_duplication_evidence_source_output = open(working_dir + '/evidences/dup_tan_source.bed', 'w')
    tandem_duplication_evidence_dest_output = open(working_dir + '/evidences/dup_tan_dest.bed', 'w')
    translocation_evidence_output = open(working_dir + '/evidences/trans.bed', 'w')
    insertion_from_evidence_output = open(working_dir + '/evidences/ins_dup.bed', 'w')

    for cluster in deletion_evidence_clusters:
        print(cluster.get_bed_entry(), file=deletion_evidence_output)
    for cluster in insertion_evidence_clusters:
        print(cluster.get_bed_entry(), file=insertion_evidence_output)
    for cluster in inversion_evidence_clusters:
        print(cluster.get_bed_entry(), file=inversion_evidence_output)
    for cluster in tandem_duplication_evidence_clusters:
        bed_entries = cluster.get_bed_entries()
        print(bed_entries[0], file=tandem_duplication_evidence_source_output)
        print(bed_entries[1], file=tandem_duplication_evidence_dest_output)
    for cluster in insertion_from_evidence_clusters:
        bed_entries = cluster.get_bed_entries()
        print(bed_entries[0], file=insertion_from_evidence_output)
        print(bed_entries[1], file=insertion_from_evidence_output)
    for translocation in completed_translocations:
        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(translocation.contig1, translocation.pos1, translocation.pos1+1, ">{0}:{1}".format(translocation.contig2, translocation.pos2), translocation.evidence, translocation.read), file=translocation_evidence_output)


def write_evidence_clusters_vcf(working_dir, clusters, genome):
    """Write evidence clusters into working directory in BED format."""
    deletion_evidence_clusters, insertion_evidence_clusters, inversion_evidence_clusters, tandem_duplication_evidence_clusters, insertion_from_evidence_clusters, completed_translocations = clusters

    if not os.path.exists(working_dir + '/evidences'):
        os.mkdir(working_dir + '/evidences')
    vcf_output = open(working_dir + '/evidences/all.vcf', 'w')

    # Write header lines
    print("##fileformat=VCFv4.3", file=vcf_output)
    print("##source=callPacV1.0", file=vcf_output)
    print("##reference={0}".format(genome), file=vcf_output)
    print("##ALT=<ID=DEL,Description=\"Deletion\">", file=vcf_output)
    print("##ALT=<ID=INV,Description=\"Inversion\">", file=vcf_output)
    print("##ALT=<ID=DUP,Description=\"Duplication\">", file=vcf_output)
    print("##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">", file=vcf_output)
    print("##ALT=<ID=INS,Description=\"Insertion\">", file=vcf_output)
    print("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">", file=vcf_output)
    print("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">", file=vcf_output)
    print("##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">", file=vcf_output)
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", file=vcf_output)

    vcf_entries = []
    for cluster in deletion_evidence_clusters:
        vcf_entries.append((cluster.get_source(), cluster.get_vcf_entry()))
    for cluster in insertion_evidence_clusters:
        vcf_entries.append((cluster.get_source(), cluster.get_vcf_entry()))
    for cluster in inversion_evidence_clusters:
        vcf_entries.append((cluster.get_source(), cluster.get_vcf_entry()))
    for cluster in tandem_duplication_evidence_clusters:
        vcf_entries.append((cluster.get_source(), cluster.get_vcf_entry()))

    # Sort and write entries to VCF
    for source, entry in sorted(vcf_entries, key=lambda pair: pair[0]):
        print(entry, file=vcf_output)

    vcf_output.close()

def write_candidates(working_dir, candidates):
    insertion_candidates, int_duplication_candidates = candidates

    if not os.path.exists(working_dir + '/candidates'):
        os.mkdir(working_dir + '/candidates')
    #deletion_candidate_output = open(working_dir + '/candidates/candidates_deletions.bed', 'w')
    insertion_candidate_source_output = open(working_dir + '/candidates/candidates_insertions_source.bed', 'w')
    insertion_candidate_dest_output = open(working_dir + '/candidates/candidates_insertions_dest.bed', 'w')
    # inversion_candidate_output = open(working_dir + '/candidates/inv.bed', 'w')
    # tandem_duplication_candidate_output = open(working_dir + '/candidates/dup_tan.bed', 'w')
    interspersed_duplication_candidate_source_output = open(working_dir + '/candidates/candidates_int_duplications_source.bed', 'w')
    interspersed_duplication_candidate_dest_output = open(working_dir + '/candidates/candidates_int_duplications_dest.bed', 'w')

    for candidate in insertion_candidates:
        bed_entries = candidate.get_bed_entries()
        print(bed_entries[0], file=insertion_candidate_source_output)
        print(bed_entries[1], file=insertion_candidate_dest_output)
    for candidate in int_duplication_candidates:
        bed_entries = candidate.get_bed_entries()
        print(bed_entries[0], file=interspersed_duplication_candidate_source_output)
        print(bed_entries[1], file=interspersed_duplication_candidate_dest_output)


def post_processing(sv_evidences, working_dir, genome):
    # Cluster evidences
    clusters = cluster_sv_evidences(sv_evidences)

    # Write evidences
    write_evidence_clusters_bed(working_dir, clusters)
    write_evidence_clusters_vcf(working_dir, clusters, genome)

    deletion_evidence_clusters, insertion_evidence_clusters, inversion_evidence_clusters, tandem_duplication_evidence_clusters, insertion_from_evidence_clusters, completed_translocations = clusters

    insertion_candidates = []
    int_duplication_candidates = []

    # Merge translocation breakpoints
    new_insertion_candidates = merge_translocations_at_deletions(completed_translocations, deletion_evidence_clusters)
    insertion_candidates.extend(new_insertion_candidates)

    new_insertion_candidates, new_int_duplication_candidates = merge_translocations_at_insertions(completed_translocations, insertion_evidence_clusters, deletion_evidence_clusters)
    insertion_candidates.extend(new_insertion_candidates)
    int_duplication_candidates.extend(new_int_duplication_candidates)

    # Merge insertions with source
    new_insertion_candidates, new_int_duplication_candidates = merge_insertions_from(insertion_from_evidence_clusters, deletion_evidence_clusters)
    insertion_candidates.extend(new_insertion_candidates)
    int_duplication_candidates.extend(new_int_duplication_candidates)

    #Write candidates
    write_candidates(working_dir, (insertion_candidates, int_duplication_candidates))