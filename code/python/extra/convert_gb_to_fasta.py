from Bio import SeqIO
import Bio
import os
import glob
import re
import sys


def createParanoidFASTAfromFile(gbkf, ext_replace='.in.fasta', gene_prefix=None):
    """
    Extracts sequence information from a FASTA (*.fasta) or GenBank file (*.gb, *.gbk *.gbff) and writes a
    simplified FASTA output.

    - *gbkf* FASTA (*.fasta) or GenBank file (*.gb, *.gbk, *. gbff) containing CDS annotations
    - *ext_replace* [default='.in.fasta'] replace the extension (as above) with this
    - *gene_prefix* [default=None] prefix gene names with this in the fasta file

    """

    if type(gbkf) == str:
        gbkf = [gbkf]

    proteins = {}
    cntr = 0
    for fasta in gbkf:
        if fasta.endswith('.fasta'):
            if cntr == 0:
                outF = fasta.replace('.fasta', '.in.fasta')
                cntr += 1
            for seq_record in SeqIO.parse(fasta, "fasta"):
                seq_record.description = ''
                if gene_prefix is not None:
                    seq_record.id = gene_prefix + seq_record.id
                proteins[seq_record.id] = seq_record

        elif fasta.endswith('.gbk') or fasta.endswith('.gb') or fasta.endswith('.gbff'):
            if cntr == 0:
                outF = fasta.replace('.gbk', '.in.fasta').replace('.gbff', '.in.fasta').replace('.gb', '.fasta')
                cntr += 1
            GBFile = file(fasta, 'r')
            GBcds = Bio.SeqIO.InsdcIO.GenBankCdsFeatureIterator(GBFile)
            for cds in GBcds:
                if cds.seq is not None:
                    cds.id = cds.name
                    cds.description = ''
                    if gene_prefix is not None:
                        cds.id = gene_prefix + cds.id
                    proteins[cds.name] = cds
            GBFile.close()
        else:
            raise RuntimeError, 'ERROR: Unknown file: {}'.format(fasta)

    print('\nProteins: {}\n'.format(len(proteins)))

    writeFASTA(outF, proteins)
    return outF


def writeFASTA(fname, sequences, paranoid_style=True):

    F = open(fname, 'w')
    cntr = 0
    #print(len(sequences))
    for seq in sequences:
        if paranoid_style:
            sequences[seq].description = ''

        F.write(sequences[seq].format('fasta'))
        cntr += 1
        if cntr >= 500:
            F.flush()
            cntr = 0
    F.flush()
    F.close()
    print('FASTA sequence file created: {}'.format(fname))

# add path of gb files
gbf = glob.glob("../../path/*.gb")

for fn in gbf:
    createParanoidFASTAfromFile(fn)

sys.exit(0)
for fn in gbf:

    input_handle = open(fn, "r")

    directory = fn[:-4]

    if not os.path.exists(directory):
        os.makedirs(directory)

    fn_out = fn.rsplit('/', 1)[-1]
    output_handle = open(directory + '/' + fn_out[:-4] + ".fasta", "w")

    for seq_record in SeqIO.parse(input_handle, "genbank"):
        print "Dealing with GenBank record %s" % seq_record.id
        for seq_feature in seq_record.features:
            if seq_feature.type == "CDS":
                assert len(seq_feature.qualifiers['translation']) == 1
                output_handle.write(">%s from %s\n%s\n" % (
                       seq_feature.qualifiers['locus_tag'][0],
                       seq_record.name,
                       seq_feature.qualifiers['translation'][0]))

    output_handle.close()
    input_handle.close()
