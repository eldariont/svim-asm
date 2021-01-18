class Candidate:
    """Candidate class for structural variant candidates. Candidates reflect the final SV types and can be merged from signatures of several reads.
    """
    def __init__(self, source_contig, source_start, source_end, genotype = "1/1"):
        self.source_contig = source_contig
        self.source_start = source_start
        self.source_end = source_end

        self.type = None
        self.genotype = genotype


    def get_source(self):
        return (self.source_contig, self.source_start, self.source_end)


    def get_key(self):
        contig, start, end = self.get_source()
        return (self.type, contig, (start + end) // 2)


    def position_distance_to(self, candidate2):
        """Return position distance between two candidates."""
        this_contig, this_start, this_end = self.get_source()
        other_contig, other_start, other_end = candidate2.get_source()
        this_center = (this_start + this_end) // 2
        other_center = (other_start + other_end) // 2
        if self.type == candidate2.type and this_contig == other_contig:
            return min(abs(this_start - other_start), abs(this_end - other_end), abs(this_center - other_center))
        else:
            return float("inf")


    def get_vcf_entry(self):
        raise NotImplementedError


class CandidateDeletion(Candidate):
    def __init__(self, source_contig, source_start, source_end, reads, bam, genotype = "1/1"):
        assert source_end >= source_start, "Deletion end ({0}:{1}) is smaller than its start ({0}:{2}). From read {3}".format(source_contig, source_end, source_start, reads)
        self.source_contig = source_contig
        contig_length = bam.get_reference_length(source_contig)
        #0-based start of the deletion (first deleted base)
        self.source_start = max(0, source_start)
        #0-based end of the deletion (one past the last deleted base)
        self.source_end = min(contig_length, source_end)

        self.type = "DEL"
        self.reads = reads
        self.genotype = genotype


    def get_vcf_entry(self, sequence_alleles = False, reference = None, read_names = False):
        contig, start, end = self.get_source()
        filters = []
        if sequence_alleles:
            ref_allele = reference.fetch(contig, max(0, start-1), end).upper()
            alt_allele = reference.fetch(contig, max(0, start-1), start).upper()
        else:
            ref_allele = "N"
            alt_allele = "<" + self.type + ">"
        info_template="SVTYPE={0};END={1};SVLEN={2}"
        info_string = info_template.format(self.type, 
                                           end, 
                                           start - end)
        if read_names:
            info_string += ";READS={0}".format(",".join(self.reads))
        return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
                    chrom=contig,
                    pos=max(1, start),
                    id="PLACEHOLDERFORID",
                    ref=ref_allele,
                    alt=alt_allele,
                    qual=".",
                    filter="PASS" if len(filters) == 0 else ";".join(filters),
                    info=info_string,
                    format="GT",
                    samples="{gt}".format(gt=self.genotype))


class CandidateInversion(Candidate):
    def __init__(self, source_contig, source_start, source_end, reads, complete, bam, genotype = "1/1"):
        assert source_end >= source_start, "Inversion end ({0}:{1}) is smaller than its start ({0}:{2}). From read {3}".format(source_contig, source_end, source_start, reads)
        self.source_contig = source_contig
        contig_length = bam.get_reference_length(source_contig)
        #0-based start of the inversion (first inverted base)
        self.source_start = max(0, source_start)
        #0-based end of the inversion (one past the last inverted base)
        self.source_end = min(contig_length, source_end)

        self.type = "INV"
        self.reads = reads
        self.complete = complete
        self.genotype = genotype

        self.complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


    def get_vcf_entry(self, sequence_alleles = False, reference = None, read_names = False):
        contig, start, end = self.get_source()
        filters = []
        if not self.complete:
            filters.append("incomplete_inversion")
        if sequence_alleles:
            ref_allele = reference.fetch(contig, start, end).upper()
            alt_allele = "".join(self.complement.get(base.upper(), base.upper()) for base in reversed(ref_allele))
        else:
            ref_allele = "N"
            alt_allele = "<" + self.type + ">"
        info_template="SVTYPE={0};END={1}"
        info_string = info_template.format(self.type, 
                                            end)
        if read_names:
            info_string += ";READS={0}".format(",".join(self.reads))
        return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
                    chrom=contig,
                    pos=start+1,
                    id="PLACEHOLDERFORID",
                    ref=ref_allele,
                    alt=alt_allele,
                    qual=".",
                    filter="PASS" if len(filters) == 0 else ";".join(filters),
                    info=info_string,
                    format="GT",
                    samples="{gt}".format(gt=self.genotype))


class CandidateInsertion(Candidate):
    def __init__(self, dest_contig, dest_start, dest_end, reads, sequence, bam, genotype = "1/1"):
        assert dest_end >= dest_start, "Insertion end ({0}:{1}) is smaller than its start ({0}:{2}). From read {3}".format(dest_contig, dest_end, dest_start, reads)
        self.dest_contig = dest_contig
        contig_length = bam.get_reference_length(dest_contig)
        #0-based start of the insertion (base after the insertion)
        self.dest_start = max(0, dest_start)
        #0-based start of the insertion (base after the insertion) + length of the insertion
        self.dest_end = min(contig_length, dest_end)

        self.type = "INS"
        self.reads = reads
        self.sequence = sequence
        self.genotype = genotype

    def get_destination(self):
        return (self.dest_contig, self.dest_start, self.dest_end)


    def get_key(self):
        return (self.type, self.dest_contig, self.dest_start)


    def get_vcf_entry(self, sequence_alleles = False, reference = None, read_names = False):
        contig, start, end = self.get_destination()
        filters = []
        if sequence_alleles:
            ref_allele = reference.fetch(contig, max(0, start-1), start).upper()
            alt_allele = ref_allele + self.sequence
        else:
            ref_allele = "N"
            alt_allele = "<" + self.type + ">"
        info_template="SVTYPE={0};END={1};SVLEN={2}"
        info_string = info_template.format(self.type, 
                                           start, 
                                           end - start) 
        if read_names:
            info_string += ";READS={0}".format(",".join(self.reads))
        return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
                    chrom=contig,
                    pos=max(1, start),
                    id="PLACEHOLDERFORID",
                    ref=ref_allele,
                    alt=alt_allele,
                    qual=".",
                    filter="PASS" if len(filters) == 0 else ";".join(filters),
                    info=info_string,
                    format="GT",
                    samples="{gt}".format(gt=self.genotype))


class CandidateDuplicationTandem(Candidate):
    def __init__(self, source_contig, source_start, source_end, copies, fully_covered, reads, bam, genotype = "1/1"):
        assert source_end >= source_start, "Tandem duplication end ({0}:{1}) is smaller than its start ({0}:{2}). From read {3}".format(source_contig, source_end, source_start, reads)
        self.source_contig = source_contig
        contig_length = bam.get_reference_length(source_contig)
        #0-based start of the region (first copied base)
        self.source_start = max(0, source_start)
        #0-based end of the region (one past the last copied base)
        self.source_end = min(contig_length, source_end)
        
        #number of additional copies
        self.copies = copies

        self.type = "DUP_TAN"
        self.reads = reads
        self.fully_covered = fully_covered
        self.genotype = genotype


    def get_destination(self):
        source_contig, source_start, source_end = self.get_source()
        return (source_contig, source_end, source_end + self.copies * (source_end - source_start))


    def get_vcf_entry_as_ins(self, sequence_alleles = False, reference = None, read_names = False):
        contig = self.source_contig
        start = self.source_start
        end = self.source_end
        svtype = "INS"
        filters = []
        if sequence_alleles:
            ref_allele = reference.fetch(contig, self.source_start, self.source_end).upper()
            alt_allele = ref_allele * (self.copies + 1)
        else:
            ref_allele = "N"
            alt_allele = "<" + svtype + ">"
        if not(self.fully_covered):
            filters.append("not_fully_covered")
        info_template="SVTYPE={0};END={1};SVLEN={2}"
        info_string = info_template.format(svtype, 
                                           end, 
                                           (end - start) * self.copies)
        if read_names:
            info_string += ";READS={0}".format(",".join(self.reads))
        return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
                    chrom=contig,
                    pos=start+1,
                    id="PLACEHOLDERFORID",
                    ref=ref_allele,
                    alt=alt_allele,
                    qual=".",
                    filter="PASS" if len(filters) == 0 else ";".join(filters),
                    info=info_string,
                    format="GT",
                    samples="{gt}".format(gt=self.genotype))


    def get_vcf_entry_as_dup(self, read_names = False):
        contig = self.source_contig
        start = self.source_start
        end = self.source_end
        length = self.source_end - self.source_start
        svtype = "DUP:TANDEM"
        filters = []
        if not(self.fully_covered):
            filters.append("not_fully_covered")
        info_template="SVTYPE={0};END={1};SVLEN={2}"
        info_string = info_template.format(svtype, 
                                           end, 
                                           length)
        if read_names:
            info_string += ";READS={0}".format(",".join(self.reads))
        return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
                    chrom=contig,
                    pos=start+1,
                    id="PLACEHOLDERFORID",
                    ref="N",
                    alt="<" + svtype + ">",
                    qual=".",
                    filter="PASS" if len(filters) == 0 else ";".join(filters),
                    info=info_string,
                    format="GT:CN",
                    samples="{gt}:{cn}".format(gt=self.genotype, cn=self.copies + 1))


class CandidateDuplicationInterspersed(Candidate):
    def __init__(self, source_contig, source_start, source_end, dest_contig, dest_start, dest_end, reads, bam, cutpaste=False, genotype = "1/1"):
        assert source_end >= source_start, "Interspersed duplication source end ({0}:{1}) is smaller than its start ({0}:{2}). From read {3}".format(source_contig, source_end, source_start, reads)
        assert dest_end >= dest_start, "Interspersed duplication destination end ({0}:{1}) is smaller than its start ({0}:{2}). From read {3}".format(dest_contig, dest_end, dest_start, reads)
        self.source_contig = source_contig
        source_contig_length = bam.get_reference_length(source_contig)
        #0-based start of the region (first copied base)
        self.source_start = max(0, source_start)
        #0-based end of the region (one past the last copied base)
        self.source_end = min(source_contig_length, source_end)

        self.dest_contig = dest_contig
        dest_contig_length = bam.get_reference_length(dest_contig)
        #0-based start of the insertion (base after the insertion)
        self.dest_start = max(0, dest_start)
        #0-based end of the insertion (base after the insertion) + length of the insertion
        self.dest_end = min(dest_contig_length, dest_end)

        self.cutpaste= cutpaste
        self.type = "DUP_INT"
        self.reads = reads
        self.genotype = genotype


    def get_destination(self):
        return (self.dest_contig, self.dest_start, self.dest_end)


    def get_key(self):
        return (self.type, self.dest_contig, self.dest_start)


    def get_vcf_entry_as_ins(self, sequence_alleles = False, reference = None, read_names = False):
        contig, start, end = self.get_destination()
        svtype = "INS"
        filters = []
        if sequence_alleles:
            ref_allele = reference.fetch(contig, max(0, start-1), start).upper()
            alt_allele = ref_allele + reference.fetch(self.source_contig, self.source_start, self.source_end).upper()
        else:
            ref_allele = "N"
            alt_allele = "<" + svtype + ">"
        info_template="SVTYPE={0};{1}END={2};SVLEN={3}"
        info_string = info_template.format(svtype, 
                                           "CUTPASTE;" if self.cutpaste else "", 
                                           start, 
                                           end - start)
        if read_names:
            info_string += ";READS={0}".format(",".join(self.reads))
        return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
                    chrom=contig,
                    pos=max(1, start),
                    id="PLACEHOLDERFORID",
                    ref=ref_allele,
                    alt=alt_allele,
                    qual=".",
                    filter="PASS" if len(filters) == 0 else ";".join(filters),
                    info=info_string,
                    format="GT",
                    samples="{gt}".format(gt=self.genotype))


    def get_vcf_entry_as_dup(self, read_names = False):
        contig, start, end = self.get_source()
        svtype = "DUP:INT"
        filters = []
        info_template="SVTYPE={0};{1}END={2};SVLEN={3}"
        info_string = info_template.format(svtype, 
                                           "CUTPASTE;" if self.cutpaste else "", 
                                           end, 
                                           end - start)
        if read_names:
            info_string += ";READS={0}".format(",".join(self.reads))
        return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
                    chrom=contig,
                    pos=start+1,
                    id="PLACEHOLDERFORID",
                    ref="N",
                    alt="<" + svtype + ">",
                    qual=".",
                    filter="PASS" if len(filters) == 0 else ";".join(filters),
                    info=info_string,
                    format="GT",
                    samples="{gt}".format(gt=self.genotype))


class CandidateBreakend(Candidate):
    def __init__(self, source_contig, source_start, source_direction, dest_contig, dest_start, dest_direction, reads, bam, genotype = "1/1"):
        if source_contig < dest_contig or (source_contig == dest_contig and source_start < dest_start):
            self.source_contig = source_contig
            source_contig_length = bam.get_reference_length(source_contig)
            #0-based source of the translocation (first base before the translocation)
            self.source_start = min(source_contig_length, max(0, source_start))
            self.source_direction = source_direction
            self.dest_contig = dest_contig
            dest_contig_length = bam.get_reference_length(dest_contig)
            #0-based destination of the translocation (first base after the translocation)
            self.dest_start = min(dest_contig_length, max(0, dest_start))
            self.dest_direction = dest_direction
        else:
            self.source_contig = dest_contig
            source_contig_length = bam.get_reference_length(dest_contig)
            #0-based source of the translocation (first base before the translocation)
            self.source_start = min(source_contig_length, max(0, dest_start))
            self.source_direction = 'fwd' if dest_direction == 'rev' else 'rev'
            self.dest_contig = source_contig
            dest_contig_length = bam.get_reference_length(source_contig)
            #0-based destination of the translocation (first base after the translocation)
            self.dest_start = min(dest_contig_length, max(0, source_start))
            self.dest_direction = 'fwd' if source_direction == 'rev' else 'rev'
        self.type = "BND"
        self.reads = reads
        self.genotype = genotype


    def get_source(self):
        return (self.source_contig, self.source_start)


    def get_destination(self):
        return (self.dest_contig, self.dest_start)
    
    def get_key(self):
        return (self.type, self.source_contig, self.source_start)

    def get_vcf_entry(self, read_names = False):
        source_contig, source_start = self.get_source()
        dest_contig, dest_start = self.get_destination()
        if (self.source_direction == 'fwd') and (self.dest_direction == 'fwd'):
            alt_string = "N[{contig}:{start}[".format(contig = dest_contig, start = dest_start+1)
        elif (self.source_direction == 'fwd') and (self.dest_direction == 'rev'):
            alt_string = "N]{contig}:{start}]".format(contig = dest_contig, start = dest_start+1)
        elif (self.source_direction == 'rev') and (self.dest_direction == 'rev'):
            alt_string = "]{contig}:{start}]N".format(contig = dest_contig, start = dest_start+1)
        elif (self.source_direction == 'rev') and (self.dest_direction == 'fwd'):
            alt_string = "[{contig}:{start}[N".format(contig = dest_contig, start = dest_start+1)
        filters = []
        info_template="SVTYPE={0}"
        info_string = info_template.format(self.type)
        if read_names:
            info_string += ";READS={0}".format(",".join(self.reads))
        return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
                    chrom=source_contig,
                    pos=source_start+1,
                    id="PLACEHOLDERFORID",
                    ref="N",
                    alt=alt_string,
                    qual=".",
                    filter="PASS" if len(filters) == 0 else ";".join(filters),
                    info=info_string,
                    format="GT",
                    samples="{gt}".format(gt=self.genotype))

    def get_vcf_entry_reverse(self, read_names = False):
        source_contig, source_start = self.get_destination()
        dest_contig, dest_start = self.get_source()
        if (self.source_direction == 'rev') and (self.dest_direction == 'rev'):
            alt_string = "N[{contig}:{start}[".format(contig = dest_contig, start = dest_start+1)
        elif (self.source_direction == 'fwd') and (self.dest_direction == 'rev'):
            alt_string = "N]{contig}:{start}]".format(contig = dest_contig, start = dest_start+1)
        elif (self.source_direction == 'fwd') and (self.dest_direction == 'fwd'):
            alt_string = "]{contig}:{start}]N".format(contig = dest_contig, start = dest_start+1)
        elif (self.source_direction == 'rev') and (self.dest_direction == 'fwd'):
            alt_string = "[{contig}:{start}[N".format(contig = dest_contig, start = dest_start+1)
        filters = []
        info_template="SVTYPE={0}"
        info_string = info_template.format(self.type)
        if read_names:
            info_string += ";READS={0}".format(",".join(self.reads))
        return "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{samples}".format(
                    chrom=source_contig,
                    pos=source_start+1,
                    id="PLACEHOLDERFORID",
                    ref="N",
                    alt=alt_string,
                    qual=".",
                    filter="PASS" if len(filters) == 0 else ";".join(filters),
                    info=info_string,
                    format="GT",
                    samples="{gt}".format(gt=self.genotype))
