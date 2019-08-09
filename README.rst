SVIM-asm - Structural variant identification method (Assembly edition)
======================================================================

SVIM-asm (pronounced *SWIM-assem*) is a structural variant caller for genome-genome alignments.
It analyzes a given sorted SAM/BAM file (preferably from minimap2) and detects five different variant classes between the query assembly and the reference: deletions, insertions, tandem and interspersed duplications and inversions.

Background
----------

.. image:: https://raw.githubusercontent.com/eldariont/svim/master/docs/SVclasses.png
    :align: center

Structural variants (SVs) are typically defined as genomic variants larger than 50bps (e.g. deletions, duplications, inversions).
Studies have shown that they affect more bases in any given genome than SNPs and small Indels taken together.
Consequently, they have a large impact on genes and regulatory regions.
This is reflected in the large number of genetic diseases that are caused by SVs.

Nowadays, SVs are usually detected using data from second-generation sequencing (Illumina) or third-generation sequencing (PacBio and Oxford Nanopore).
Typically, the reads from a sequencing experiment are first aligned to a reference genome before the alignments are analyzed for characteristic signatures of SVs.
Recently, substantial advances in sequencing technology and software development have made the de novo assembly of large mammalian genomes more efficient than ever.
Accurate assemblies of the human genome can now be generated in a few days and at a fraction of its former cost. `[Shafin et al.] <https://www.biorxiv.org/content/10.1101/715722v1>`_

Similarly to raw sequencing reads, the genome assemblies can be aligned to another genome to uncover genomic rearrangements and structural variants.
Our tool, SVIM-asm, detects structural variants between different assemblies or reference genomes from given genome-genome alignments.
It is fast (<5 min for a human genome-genome alignment), easy to use and detects all major variant types.

Installation
------------

.. code-block:: bash

    #Install from github (requires Python 3)
    git clone https://github.com/eldariont/svim-asm.git
    cd svim-asm
    pip3 install .

Execution
-----

SVIM-asm analyzes alignments between a query assembly and a reference assembly in SAM/BAM format. 
We recommend to produce the alignments using `minimap2 <https://github.com/lh3/minimap2>`_:

.. code-block:: bash

    minimap2 --paf-no-hit -a -x asm5 --cs -r2k -t <num_threads> <reference.fa> <assembly.fasta> > <alignments.sam>
    samtools sort -m4G -@4 -o <alignments.sorted.bam> <alignments.sam>
    svim-asm alignment <working_dir> <alignments.sorted.bam> <reference.fa>

Output
------

SVIM-asm creates all output files in the given working directory.
The following files are produced:

- ``variants.vcf`` contains the detected SVs in VCF format (see http://samtools.github.io/hts-specs/VCFv4.2.pdf)
- ``sv-lengths.png`` contains a histogram of SV sizes
- ``SVIM_<day>_<time>.log`` contains the same logging output as the command line 

Contact
-------

If you experience problems or have suggestions please create an issue or a pull request or contact heller_d@molgen.mpg.de.

Citation
---------

SVIM-asm is a fork of our long-read caller SVIM. Feel free to read and cite our paper in Bioinformatics: https://doi.org/10.1093/bioinformatics/btz041

License
-------

The project is licensed under the GNU General Public License.
