import sys
import os
import logging
import argparse


def parse_arguments(program_version, arguments = sys.argv[1:]):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""SVIM-asm (pronounced SWIM-assem) is a structural variant caller for genome-genome alignments. 
It discriminates five different variant classes: deletions, insertions, tandem and interspersed duplications and inversions.
SVIM-asm analyzes alignments between a haploid or diploid query assembly and a reference assembly in SAM/BAM format. 
We recommend to produce the alignments using minimap2.

SVIM-asm has an haploid and a diploid mode depending on the input assembly and performs the following steps:
- COLLECT detects SVs from genome-genome alignments in BAM format
- PAIR merges the SV calls from the two haplotypes of a diploid assembly (diploid mode only)
- OUTPUT prints the found SVs in VCF format
""")

    subparsers = parser.add_subparsers(help='modes', dest='sub')
    parser.add_argument('--version',
                        '-v',
                        action='version',
                        version='%(prog)s {version}'.format(version=program_version))
    
    parser_haploid = subparsers.add_parser('haploid',
                                        help='Detect SVs from the alignment of an haploid query assembly to a reference assembly')
    parser_haploid.add_argument('working_dir',
                             type=os.path.abspath,
                             help='Working and output directory. \
                                   Existing files in the directory are overwritten. \
                                   If the directory does not exist, it is created.')
    parser_haploid.add_argument('bam_file',
                             type=str,
                             help='SAM/BAM file with alignment of query assembly to reference assembly (needs to be coordinate-sorted and indexed)')
    parser_haploid.add_argument('genome',
                               type=str,
                               help='Reference genome file that the assembly was aligned to (FASTA)')
    parser_haploid.add_argument('--verbose',
                              action='store_true',
                              help='Enable more verbose logging (default: %(default)s)')
    group_haploid_collect = parser_haploid.add_argument_group('COLLECT')
    group_haploid_collect.add_argument('--min_mapq',
                                      type=int,
                                      default=20,
                                      help='Minimum mapping quality of alignments to consider (default: %(default)s). \
                                            Alignments with a lower mapping quality are ignored.')
    group_haploid_collect.add_argument('--min_sv_size',
                                      type=int,
                                      default=40,
                                      help='Minimum SV size to detect (default: %(default)s). \
                                            SVIM can potentially detect events of any size but is limited by the \
                                            signal-to-noise ratio in the input alignments. That means that more \
                                            accurate assemblies and alignments enable the detection of smaller events.')
    group_haploid_collect.add_argument('--max_sv_size',
                                      type=int,
                                      default=100000,
                                      help='Maximum SV size to detect (default: %(default)s). \
                                              This parameter is used to distinguish long deletions (and inversions) from \
                                              translocations which cannot be distinguished from the alignment alone. \
                                              Split read segments mapping far apart on the reference could either \
                                              indicate a very long deletion (inversion) or a translocation breakpoint. \
                                              SVIM calls a translocation breakpoint if the mapping distance is larger \
                                              than this parameter and a deletion (or inversion) if it is smaller or equal.')
    group_haploid_collect.add_argument('--query_gap_tolerance',
                                      type=int,
                                      default=50,
                                      help='Maximum tolerated gap between adjacent alignment segments on the query \
                                            (default: %(default)s). Example: \
                                            Deletions are detected from two subsequent segments of a split query sequence that are mapped \
                                            far apart from each other on the reference. The query gap tolerance determines \
                                            the maximum tolerated length of the query gap between both segments. If there is an \
                                            unaligned query segment larger than this value between the two segments, no deletion is called.')
    group_haploid_collect.add_argument('--query_overlap_tolerance',
                                      type=int,
                                      default=50,
                                      help='Maximum tolerated overlap between adjacent alignment segments on the query \
                                            (default: %(default)s). Example: \
                                            Deletions are detected from two subsequent segments of a split query sequence that are mapped \
                                            far apart from each other on the reference. The query overlap tolerance determines \
                                            the maximum tolerated length of an overlap between both segments in the query. If the \
                                            overlap between the two segments in the query is larger than this value, no deletion is called.')
    group_haploid_collect.add_argument('--reference_gap_tolerance',
                                      type=int,
                                      default=50,
                                      help='Maximum tolerated gap between adjacent alignment segments on the reference \
                                            (default: %(default)s). Example: \
                                            Insertions are detected from two segments of a split query sequence that are mapped \
                                            right next to each other on the reference but with unaligned sequence between them on the query. \
                                            The reference gap tolerance determines the maximum tolerated length of the reference gap between \
                                            both segments. If there is a reference gap larger than this value between the two segments, no \
                                            insertion is called.')
    group_haploid_collect.add_argument('--reference_overlap_tolerance',
                                      type=int,
                                      default=50,
                                      help='Maximum tolerated overlap between adjacent alignment segments on the reference \
                                            (default: %(default)s). Example: \
                                            Insertions are detected from two segments of a split query sequence that are mapped \
                                            right next to each other on the reference but with unaligned sequence between them on the query. \
                                            The reference overlap tolerance determines the maximum tolerated length of an overlap between \
                                            both segments on the reference. If there is a reference gap larger than this value between the \
                                            two segments, no insertion is called.')
    
    group_haploid_output = parser_haploid.add_argument_group('OUTPUT')
    group_haploid_output.add_argument('--sample',
                                        type=str,
                                        default="Sample",
                                        help='Sample ID to include in output vcf file (default: %(default)s)')
    group_haploid_output.add_argument('--types',
                                        type=str,
                                        default="DEL,INS,INV,DUP:TANDEM,DUP:INT,BND",
                                        help='SV types to include in output VCF (default: %(default)s). \
                                              Give a comma-separated list of SV types. The possible SV types are: DEL (deletions), \
                                              INS (novel insertions), INV (inversions), DUP:TANDEM (tandem duplications), \
                                              DUP:INT (interspersed duplications), BND (breakends).')
    group_haploid_output.add_argument('--symbolic_alleles',
                                        action='store_true',
                                        help='Use symbolic alleles, such as <DEL> or <INV> in the VCF output (default: %(default)s). \
                                              By default, deletions, insertions, and inversions are represented by their nucleotide sequence in the output VCF.')
    group_haploid_output.add_argument('--tandem_duplications_as_insertions',
                                        action='store_true',
                                        help='Represent tandem duplications as insertions in output VCF (default: %(default)s). \
                                              By default, tandem duplications are represented by the SVTYPE=DUP:TANDEM and the genomic source is given by the \
                                              POS and END tags. When enabling this option, duplications are instead represented by the SVTYPE=INS \
                                              and POS and END both give the insertion point of the duplication.')
    group_haploid_output.add_argument('--interspersed_duplications_as_insertions',
                                        action='store_true',
                                        help='Represent interspersed duplications as insertions in output VCF (default: %(default)s). \
                                              By default, interspersed duplications are represented by the SVTYPE=DUP:INT and the genomic source is given by the \
                                              POS and END tags. When enabling this option, duplications are instead represented by the SVTYPE=INS \
                                              and POS and END both give the insertion point of the duplication.')
    group_haploid_output.add_argument('--query_names',
                                        action='store_true',
                                        help='Output names of supporting query sequences in INFO tag of VCF (default: %(default)s). \
                                              If enabled, the INFO/READS tag contains the list of names of the supporting query sequences.')

    parser_diploid = subparsers.add_parser('diploid',
                                        help='Detect SVs from the alignment of a diploid query assembly to a reference assembly')
    parser_diploid.add_argument('working_dir',
                             type=os.path.abspath,
                             help='Working and output directory. \
                                   Existing files in the directory are overwritten. \
                                   If the directory does not exist, it is created.')
    parser_diploid.add_argument('bam_file1',
                             type=str,
                             help='SAM/BAM file with alignment of query assembly\'s first haplotype to reference assembly (needs to be coordinate-sorted and indexed)')
    parser_diploid.add_argument('bam_file2',
                             type=str,
                             help='SAM/BAM file with alignment of query assembly\'s second haplotype to reference assembly (needs to be coordinate-sorted and indexed)')
    parser_diploid.add_argument('genome',
                               type=str,
                               help='Reference genome file that the assembly was aligned to (FASTA)')
    parser_diploid.add_argument('--verbose',
                              action='store_true',
                              help='Enable more verbose logging (default: %(default)s)')
    group_diploid_collect = parser_diploid.add_argument_group('COLLECT')
    group_diploid_collect.add_argument('--min_mapq',
                                      type=int,
                                      default=20,
                                      help='Minimum mapping quality of alignments to consider (default: %(default)s). \
                                            Alignments with a lower mapping quality are ignored.')
    group_diploid_collect.add_argument('--min_sv_size',
                                      type=int,
                                      default=40,
                                      help='Minimum SV size to detect (default: %(default)s). \
                                            SVIM can potentially detect events of any size but is limited by the \
                                            signal-to-noise ratio in the input alignments. That means that more \
                                            accurate assemblies and alignments enable the detection of smaller events.')
    group_diploid_collect.add_argument('--max_sv_size',
                                      type=int,
                                      default=100000,
                                      help='Maximum SV size to detect (default: %(default)s). \
                                              This parameter is used to distinguish long deletions (and inversions) from \
                                              translocations which cannot be distinguished from the alignment alone. \
                                              Split read segments mapping far apart on the reference could either \
                                              indicate a very long deletion (inversion) or a translocation breakpoint. \
                                              SVIM calls a translocation breakpoint if the mapping distance is larger \
                                              than this parameter and a deletion (or inversion) if it is smaller or equal.')
    group_diploid_collect.add_argument('--query_gap_tolerance',
                                      type=int,
                                      default=50,
                                      help='Maximum tolerated gap between adjacent alignment segments on the query \
                                            (default: %(default)s). Example: \
                                            Deletions are detected from two subsequent segments of a split query sequence that are mapped \
                                            far apart from each other on the reference. The query gap tolerance determines \
                                            the maximum tolerated length of the query gap between both segments. If there is an \
                                            unaligned query segment larger than this value between the two segments, no deletion is called.')
    group_diploid_collect.add_argument('--query_overlap_tolerance',
                                      type=int,
                                      default=50,
                                      help='Maximum tolerated overlap between adjacent alignment segments on the query \
                                            (default: %(default)s). Example: \
                                            Deletions are detected from two subsequent segments of a split query sequence that are mapped \
                                            far apart from each other on the reference. The query overlap tolerance determines \
                                            the maximum tolerated length of an overlap between both segments in the query. If the \
                                            overlap between the two segments in the query is larger than this value, no deletion is called.')
    group_diploid_collect.add_argument('--reference_gap_tolerance',
                                      type=int,
                                      default=50,
                                      help='Maximum tolerated gap between adjacent alignment segments on the reference \
                                            (default: %(default)s). Example: \
                                            Insertions are detected from two segments of a split query sequence that are mapped \
                                            right next to each other on the reference but with unaligned sequence between them on the query. \
                                            The reference gap tolerance determines the maximum tolerated length of the reference gap between \
                                            both segments. If there is a reference gap larger than this value between the two segments, no \
                                            insertion is called.')
    group_diploid_collect.add_argument('--reference_overlap_tolerance',
                                      type=int,
                                      default=50,
                                      help='Maximum tolerated overlap between adjacent alignment segments on the reference \
                                            (default: %(default)s). Example: \
                                            Insertions are detected from two segments of a split query sequence that are mapped \
                                            right next to each other on the reference but with unaligned sequence between them on the query. \
                                            The reference overlap tolerance determines the maximum tolerated length of an overlap between \
                                            both segments on the reference. If there is a reference gap larger than this value between the \
                                            two segments, no insertion is called.')

    group_diploid_pair = parser_diploid.add_argument_group('PAIR')
    group_diploid_pair.add_argument('--partition_max_distance',
                                        type=int,
                                        default=1000,
                                        help='Maximum distance in bp between SVs in a partition (default: %(default)s). \
                                              Before pairing, the SV signatures are divided into coarse partitions. This parameter \
                                              determines the maximum distance between two subsequent signatures in the same partition. \
                                              If the distance between two subsequent signatures \
                                              is larger than this parameter, they are distributed into separate partitions.')
    group_diploid_pair.add_argument('--max_edit_distance',
                                        type=int,
                                        default=200,
                                        help='Maximum edit distance between both alleles to be paired up into a homozygous call (default: %(default)s).')

    group_diploid_output = parser_diploid.add_argument_group('OUTPUT')
    group_diploid_output.add_argument('--sample',
                                        type=str,
                                        default="Sample",
                                        help='Sample ID to include in output vcf file (default: %(default)s)')
    group_diploid_output.add_argument('--types',
                                        type=str,
                                        default="DEL,INS,INV,DUP:TANDEM,DUP:INT,BND",
                                        help='SV types to include in output VCF (default: %(default)s). \
                                              Give a comma-separated list of SV types. The possible SV types are: DEL (deletions), \
                                              INS (novel insertions), INV (inversions), DUP:TANDEM (tandem duplications), \
                                              DUP:INT (interspersed duplications), BND (breakends).')
    group_diploid_output.add_argument('--symbolic_alleles',
                                        action='store_true',
                                        help='Use symbolic alleles, such as <DEL> or <INV> in the VCF output (default: %(default)s). \
                                              By default, deletions, insertions, and inversions are represented by their nucleotide sequence in the output VCF.')
    group_diploid_output.add_argument('--tandem_duplications_as_insertions',
                                        action='store_true',
                                        help='Represent tandem duplications as insertions in output VCF (default: %(default)s). \
                                              By default, tandem duplications are represented by the SVTYPE=DUP:TANDEM and the genomic source is given by the \
                                              POS and END tags. When enabling this option, duplications are instead represented by the SVTYPE=INS \
                                              and POS and END both give the insertion point of the duplication.')
    group_diploid_output.add_argument('--interspersed_duplications_as_insertions',
                                        action='store_true',
                                        help='Represent interspersed duplications as insertions in output VCF (default: %(default)s). \
                                              By default, interspersed duplications are represented by the SVTYPE=DUP:INT and the genomic source is given by the \
                                              POS and END tags. When enabling this option, duplications are instead represented by the SVTYPE=INS \
                                              and POS and END both give the insertion point of the duplication.')
    group_diploid_output.add_argument('--query_names',
                                        action='store_true',
                                        help='Output names of supporting query sequences in INFO tag of VCF (default: %(default)s). \
                                              If enabled, the INFO/READS tag contains the list of names of the supporting query sequences.')
    return parser.parse_args(arguments)
