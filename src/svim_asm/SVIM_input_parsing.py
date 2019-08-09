import sys
import os
import logging
import argparse


def parse_arguments(program_version, arguments = sys.argv[1:]):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""SVIM-asm (pronounced SWIM-assem) is a structural variant caller for genome-genome alignments. 
It discriminates five different variant classes: deletions, insertions, tandem and interspersed duplications and inversions.

SVIM consists of two major steps:
- COLLECT detects SVs in genome-genome alignments in BAM format
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
                             help='SAM/BAM file with alignment of query assembly to reference assembly (needs to be coordinate-sorted)')
    parser_haploid.add_argument('genome',
                               type=str,
                               help='Reference genome file that the assembly was aligned to (FASTA)')
    group_haploid_collect = parser_haploid.add_argument_group('COLLECT')
    group_haploid_collect.add_argument('--min_mapq',
                                      type=int,
                                      default=20,
                                      help='Minimum mapping quality of reads to consider (default: %(default)s). \
                                            Reads with a lower mapping quality are ignored.')
    group_haploid_collect.add_argument('--min_sv_size',
                                      type=int,
                                      default=40,
                                      help='Minimum SV size to detect (default: %(default)s). \
                                            SVIM can potentially detect events of any size but is limited by the \
                                            signal-to-noise ratio in the input alignments. That means that more \
                                            accurate reads and alignments enable the detection of smaller events. \
                                            For current PacBio or Nanopore data, we would recommend a minimum size \
                                            of 40bp or larger.')
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
    group_haploid_collect.add_argument('--segment_gap_tolerance',
                                      type=int,
                                      default=10,
                                      help='Maximum tolerated gap between adjacent alignment segments (default: %(default)s). \
                                            This parameter applies to gaps on the reference and the read. Example: \
                                            Deletions are detected from two subsequent segments of a split read that are mapped \
                                            far apart from each other on the reference. The segment gap tolerance determines \
                                            the maximum tolerated length of the read gap between both segments. If there is an \
                                            unaligned read segment larger than this value between the two segments, no deletion is called.')
    group_haploid_collect.add_argument('--segment_overlap_tolerance',
                                      type=int,
                                      default=5,
                                      help='Maximum tolerated overlap between adjacent alignment segments (default: %(default)s). \
                                            This parameter applies to overlaps on the reference and the read. Example: \
                                            Deletions are detected from two subsequent segments of a split read that are mapped \
                                            far apart from each other on the reference. The segment overlap tolerance determines \
                                            the maximum tolerated length of an overlap between both segments on the read. If the \
                                            overlap between the two segments on the read is larger than this value, no deletion is called.')
    
    group_haploid_output = parser_haploid.add_argument_group('OUTPUT')
    group_haploid_output.add_argument('--sample',
                                        type=str,
                                        default="Sample",
                                        help='Sample ID to include in output vcf file (default: %(default)s)')
    group_haploid_output.add_argument('--types',
                                        type=str,
                                        default="DEL,INS,INV,DUP_TAN,DUP_INT,BND",
                                        help='SV types to include in output VCF (default: %(default)s). \
                                              Give a comma-separated list of SV types. The possible SV types are: DEL (deletions), \
                                              INS (novel insertions), INV (inversions), DUP_TAN (tandem duplications), \
                                              DUP_INT (interspersed duplications), BND (breakends).')
    group_haploid_output.add_argument('--symbolic_alleles',
                                        action='store_true',
                                        help='Use symbolic alleles, such as <DEL> or <INV> in the VCF output (default: %(default)s). \
                                              By default, deletions, insertions, and inversions are represented by their nucleotide sequence in the output VCF.')
    group_haploid_output.add_argument('--duplications_as_insertions',
                                        action='store_true',
                                        help='Represent tandem and interspersed duplications as insertions in output VCF (default: %(default)s). \
                                              By default, duplications are represented by the SVTYPE=DUP and the genomic source is given by the \
                                              POS and END tags. When enabling this option, duplications are instead represented by the SVTYPE=INS \
                                              and POS and END both give the insertion point of the duplication.')
    group_haploid_output.add_argument('--read_names',
                                        action='store_true',
                                        help='Output names of supporting reads in INFO tag of VCF (default: %(default)s). \
                                              If enabled, the INFO/READS tag contains the list of names of the supporting reads.')

    parser_diploid = subparsers.add_parser('diploid',
                                        help='Detect SVs from the alignment of a diploid query assembly to a reference assembly')
    parser_diploid.add_argument('working_dir',
                             type=os.path.abspath,
                             help='Working and output directory. \
                                   Existing files in the directory are overwritten. \
                                   If the directory does not exist, it is created.')
    parser_diploid.add_argument('bam_file1',
                             type=str,
                             help='SAM/BAM file with alignment of query assembly\'s first haplotype to reference assembly (needs to be coordinate-sorted)')
    parser_diploid.add_argument('bam_file2',
                             type=str,
                             help='SAM/BAM file with alignment of query assembly\'s second haplotype to reference assembly (needs to be coordinate-sorted)')
    parser_diploid.add_argument('genome',
                               type=str,
                               help='Reference genome file that the assembly was aligned to (FASTA)')
    group_diploid_collect = parser_diploid.add_argument_group('COLLECT')
    group_diploid_collect.add_argument('--min_mapq',
                                      type=int,
                                      default=20,
                                      help='Minimum mapping quality of reads to consider (default: %(default)s). \
                                            Reads with a lower mapping quality are ignored.')
    group_diploid_collect.add_argument('--min_sv_size',
                                      type=int,
                                      default=40,
                                      help='Minimum SV size to detect (default: %(default)s). \
                                            SVIM can potentially detect events of any size but is limited by the \
                                            signal-to-noise ratio in the input alignments. That means that more \
                                            accurate reads and alignments enable the detection of smaller events. \
                                            For current PacBio or Nanopore data, we would recommend a minimum size \
                                            of 40bp or larger.')
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
    group_diploid_collect.add_argument('--segment_gap_tolerance',
                                      type=int,
                                      default=10,
                                      help='Maximum tolerated gap between adjacent alignment segments (default: %(default)s). \
                                            This parameter applies to gaps on the reference and the read. Example: \
                                            Deletions are detected from two subsequent segments of a split read that are mapped \
                                            far apart from each other on the reference. The segment gap tolerance determines \
                                            the maximum tolerated length of the read gap between both segments. If there is an \
                                            unaligned read segment larger than this value between the two segments, no deletion is called.')
    group_diploid_collect.add_argument('--segment_overlap_tolerance',
                                      type=int,
                                      default=5,
                                      help='Maximum tolerated overlap between adjacent alignment segments (default: %(default)s). \
                                            This parameter applies to overlaps on the reference and the read. Example: \
                                            Deletions are detected from two subsequent segments of a split read that are mapped \
                                            far apart from each other on the reference. The segment overlap tolerance determines \
                                            the maximum tolerated length of an overlap between both segments on the read. If the \
                                            overlap between the two segments on the read is larger than this value, no deletion is called.')
    
    group_diploid_output = parser_diploid.add_argument_group('OUTPUT')
    group_diploid_output.add_argument('--sample',
                                        type=str,
                                        default="Sample",
                                        help='Sample ID to include in output vcf file (default: %(default)s)')
    group_diploid_output.add_argument('--types',
                                        type=str,
                                        default="DEL,INS,INV,DUP_TAN,DUP_INT,BND",
                                        help='SV types to include in output VCF (default: %(default)s). \
                                              Give a comma-separated list of SV types. The possible SV types are: DEL (deletions), \
                                              INS (novel insertions), INV (inversions), DUP_TAN (tandem duplications), \
                                              DUP_INT (interspersed duplications), BND (breakends).')
    group_diploid_output.add_argument('--symbolic_alleles',
                                        action='store_true',
                                        help='Use symbolic alleles, such as <DEL> or <INV> in the VCF output (default: %(default)s). \
                                              By default, deletions, insertions, and inversions are represented by their nucleotide sequence in the output VCF.')
    group_diploid_output.add_argument('--duplications_as_insertions',
                                        action='store_true',
                                        help='Represent tandem and interspersed duplications as insertions in output VCF (default: %(default)s). \
                                              By default, duplications are represented by the SVTYPE=DUP and the genomic source is given by the \
                                              POS and END tags. When enabling this option, duplications are instead represented by the SVTYPE=INS \
                                              and POS and END both give the insertion point of the duplication.')
    group_diploid_output.add_argument('--read_names',
                                        action='store_true',
                                        help='Output names of supporting reads in INFO tag of VCF (default: %(default)s). \
                                              If enabled, the INFO/READS tag contains the list of names of the supporting reads.')
    return parser.parse_args(arguments)