#
# Fix BLAST report chains from PDB to MMCIF (PDB)
#

from temp2 import pdb_proteins
import sys

def main():
    pdb_path = "../Workflow/Interfaces/"

    print "Loading PDB-Uniprot mapping..."
    pdb_uniprot, chain_author_to_mmcif = pdb_proteins(pdb_path)
    # print chain_author_to_mmcif
    # sys.exit(0)

    with open("human_deltablast_report.tab", 'r') as f, open("human_deltablast_report2.tab", 'w') as o:
        # for line in islice(f, 1, None):
        #query   hit     score   bit_score       evalue  identical       positive        gaps    align_len       qfrom   qto     hfrom   hto     qseq    hseq
        for i, line in enumerate(f):
            if i == 0:
                o.write(line)
                continue
            fields = line.split("\t")
            hit = fields[1]
            pdb, chain_author = hit.split("_")
            chain_mmcif = chain_author_to_mmcif.get((pdb, chain_author), None)
            # print chain_author, chain_mmcif
            if chain_mmcif is None:
                print "PDB chain not found in Biounit PDB-MMCIF mapping. Skipping BLAST hit", pdb, chain_author
                continue
                # chain_mmcif = chain_author
            fields = list(fields)
            fields[1] = pdb + "|" + chain_mmcif
            o.write("\t".join(fields))

if __name__ == '__main__':
    main()

    # Not so many hits skipped
    # 6816762 human_deltablast_report.tab
    # 6656019 human_deltablast_report2.tab
    #
