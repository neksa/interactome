from collections import defaultdict

import numpy as np
import textwrap

from ..structures.complexes import Complexes
from ..structures.SIFTS import SIFTS, PDBChain
from ..structures.mmcif import mmCifFile

""" Creates annotations of interactions for Jalview

Analyze interaction of  nucleosome core particles with other proteins
Create annotation per structure per chain (in PDB notation) in Jalview format

Count the number of heavy atoms interacting for in each residue in the PDB sequence
"""

PDB_MAPPING = "tests/data/mapping/pdb_chain_uniprot_201807.tsv"

NUCLEOSOMES_FILE = "tests/data/nucleosomes/ncp_pdb.csv"
UNIPROT_MAPPING = "tests/data/nucleosomes/uniprot.tab"
RESULTS_FILE = "tests/data/nucleosomes/ncp_{}.jalview"
ALIGNMENT_FILE = "tests/data/nucleosomes/ncp_{}.fasta"


def read_name_mapping():
    name_mapping = """
    RCC1-3MVDC
    RCC1-3MVDD
    BMI-1-4R8PA
    BMI-1-4R8PC
    BMI-1-4R8PD
    BMI-1-4R8PB
    Sir3-3TU4A
    Sir3-3TU4B
    Sir3-3TU4C
    Sir3-3TU4D
    Snf2-5X0YA
    Snf2-5X0YB
    GAG-5MLUA
    GAG-5MLUC
    GAG-5MLUD
    Chd1-5O9GA
    Chd1-5O9GB
    CENP-N-6C0WA
    CENP-N-6C0WB
    CENP-N-6C0WD
    CENP-C-4X23C
    CENP-C-4X23D
    CENP-C-4X23E
    CENP-C-4X23F
    CENP-C-4X23G
    CENP-C-4X23H
    CENP-C-4X23K
    CENP-C-4X23L
    INO80-6FMLT
    LANA-1ZLAG
    LANA-1ZLAH
    ACTR5-6FMLQ
    ACTR5-6FMLS
    ORF73-5GTCG
    ORF73-5GTCH
    Set8-5HQ2B
    SRM1-5HQ2H
    SGF11-4ZUXD
    SGF11-4ZUXC
    UL123-5E5AG
    UL123-5E5AH
    UBP8-4ZUXB
    Ubiquitin-4ZUXA
    TP53BP1/Ubiquitin-5KGFC
    TP53BP1-5KGFB
    TP53BP1-5KGFD
    Canonical-1KX5C
    Canonical-1KX5A
    Canonical-1KX5B
    Canonical-1KX5D
    6FQ8C-vs-6FQ5C
    6FQ8G-vs-6FQ5G
    6FQ8B-vs-6FQ5B
    6FQ8F-vs-6FQ5F
    6FQ8A-vs-6FQ5A
    6FQ8E-vs-6FQ5E
    6FQ8D-vs-6FQ5D
    6FQ8H-vs-6FQ5H
    """
    mapping = {}
    for n in name_mapping.splitlines():
        name = n.strip()
        pdb_chain = name[-5:]
        pdb = pdb_chain[-5:-1]
        chain = pdb_chain[-1:]
        mapping[(pdb.upper(), chain)] = name
    return mapping


def read_chain_mapping(fname):
    """
    pdb     chain   chain_author    uniprot seq_aln_begin   seq_aln_begin_ins       seq_aln_end     seq_aln_end_ins db_aln_begin    db_aln_end      auth_aln_begin  auth_aln_end
    200l    A       A       P00720  1       ?       164     ?       1       164     1       164
    101m    A       A       P02185  2       ?       154     ?       1       153     1       153    """

    chain_mapping = {}
    with open(fname) as f:
        for line in f:
            pdb, chain, chain_author = line.strip().split()[:3]
            chain_mapping[(pdb, chain)] = chain_author
    return chain_mapping


# def read_uniprot_mapping(fname):
#     mapping = {}
#     with open(fname) as f:
#         for line in f:
#             if line.startswith("#"):
#                 continue
#             _, uniprot, uniprot_name, _, protein_name, gene, organism = line.strip().split("\t")[:7]
#             mapping[uniprot] = protein_name
#     return mapping

def read_uniprot_mapping(fname):
    mapping = {}
    with open(fname) as f:
        for line in f:
            if line.startswith("Entry"):
                continue
            uniprot, uniprot_name, _, protein_name, gene, organism = line.strip().split("\t")[:6]
            mapping[uniprot] = protein_name
    return mapping


def main():
    # d = "/net/pan1/interactomes/pipeline/"
    d = "/panfs/pan1.be-md.ncbi.nlm.nih.gov/interactomes/pipeline/"

    # roc_d = d + "/Benchmarks/Ecoli"
    # matches_ecoli = d + "/Alignments/matches_ecoli.tab"
    # pairs = count_scored_pairs(matches_ecoli, 0.25)  # , 3.0)
    # p = set()
    # for pair in pairs:
    #     p.add(pair[0])
    #     p.add(pair[1])
    # print "N proteins = ", len(p)
    # print "N interactions = ", len(pairs)

    pdb_structures = d + "../pdb"
    pdb_interfaces = d + "Interactome/Workflow/Interfaces/{}/{}_atomic_contacts_5.0A.tab"
    pdb_templates = d + "Interactome/Workflow/Structures/pdb_templates_5A.tab"
    pdb_proteins = d + "Interactome/Workflow/Structures/pdb_proteins.tab"
    pdb_mapping_fname = d + "interactome/" + PDB_MAPPING
    uniprot_mapping_fname = d + "interactome/" + UNIPROT_MAPPING
    nucleosomes = d + "interactome/" + NUCLEOSOMES_FILE

    selected_pdbs = defaultdict(list)
    with open(nucleosomes) as f:
        for line in f:
            if len(line.strip()) > 0:
                pdb_chain = line.strip().split()[0].lower()
                pdb = pdb_chain[:4]
                chain = pdb_chain[4:]
                selected_pdbs[pdb].append(chain.upper())
    # print(len(selected_pdbs))
    # print(selected_pdbs)

    chain_mapping = read_chain_mapping(pdb_proteins)
    uniprot_mapping = read_uniprot_mapping(uniprot_mapping_fname)
    # print(uniprot_mapping)

    complexes = Complexes()
    templates = complexes.loadTemplates(pdb_templates) # , mapping)

    sifts = SIFTS()
    sifts_mapping = sifts.getPDBMapping(pdb_mapping_fname)

    name_mapping = read_name_mapping()

    results = {}
    sequences = {}

    for template in templates:
        (pdb, chainA_orig, chainB_orig), (site1, site2) = template
        chainA = chainA_orig.split("_")[0]
        chainB = chainB_orig.split("_")[0]

        if pdb in selected_pdbs.keys():
            chainA_author = chain_mapping[(pdb, chainA)]
            chainB_author = chain_mapping[(pdb, chainB)]
            if chainA_author not in selected_pdbs[pdb] and chainB_author not in selected_pdbs[pdb]:
                continue
            # print(pdb, chainA_author, chainB_author)

            chainA_uniprot = sifts_mapping.get((pdb, chainA_author), pdb + "_" + chainA_author)
            if type(chainA_uniprot) is PDBChain:
                chainA_uniprot_protein = chainA_uniprot.protein
                # print(chainA_uniprot_protein)
                chainA_uniprot = uniprot_mapping.get(chainA_uniprot_protein, chainA_uniprot_protein)

            chainB_uniprot = sifts_mapping.get((pdb, chainB_author), pdb + "_" + chainB_author)
            if type(chainB_uniprot) is PDBChain:
                chainB_uniprot_protein = chainB_uniprot.protein
                # print(chainB_uniprot_protein)
                chainB_uniprot = uniprot_mapping.get(chainB_uniprot_protein, chainB_uniprot_protein)

            # Only consider histone proteins interacting with non-histone proteins
            histone = None
            chain = None
            if chainA_uniprot.startswith("His") and not chainB_uniprot.startswith("His"):
                if "H3" in chainA_uniprot:
                    histone = "H3"
                if "H4" in chainA_uniprot:
                    histone = "H4"
                if "H2A" in chainA_uniprot:
                    histone = "H2A"
                if "H2B" in chainA_uniprot:
                    histone = "H2B"

                if histone is None:
                    print("Unidentified histone", chainA_uniprot)
                chain = chainA_author

            elif chainB_uniprot.startswith("His") and not chainA_uniprot.startswith("His"):
                if "H3" in chainB_uniprot:
                    histone = "H3"
                if "H4" in chainB_uniprot:
                    histone = "H4"
                if "H2A" in chainB_uniprot:
                    histone = "H2A"
                if "H2B" in chainB_uniprot:
                    histone = "H2B"

                if histone is None:
                    print("Unidentified histone", chainB_uniprot)
                chain = chainB_author
            else:
                pass
                # print(chainA_uniprot, chainB_uniprot)

            # print(chain, histone)

            if histone and chain:
                print(pdb, chainA_author, chainA_uniprot_protein, chainA_uniprot, chainB_author, chainB_uniprot_protein, chainB_uniprot)
                fname_int = pdb_interfaces.format(pdb[1:3], pdb)
                # print(fname_int)
                contacts = complexes.getInterface(fname_int, 5.0)
                print(contacts[(chainA_orig, chainB_orig)][0])

                s = mmCifFile(pdb_structures, pdb)
                residue_values = {}
                max_resi = 0
                min_resi = 1000
                seq = ""
                for a in s.iterAtoms():
                    # if pdb.upper() == "5HQ2":
                    #     print(a.chain_author)
                    if a.chain_author != chain:
                        continue
                    if a.resi in residue_values:
                        continue
                    max_resi = max(max_resi, a.resi)
                    min_resi = min(min_resi, a.resi)
                    seq += a.resn_short
                    # print(a.chain_author, a.resi, a.resn_short)
                    residue_values[a.resi] = contacts[(chainA_orig, chainB_orig)][0 if chainA_author == chain else 1][(a.resn_short, str(a.resi))]

                sequences[(pdb, chain)] = seq

                offset = min_resi - 1

                # print(residue_values)

                values = []
                for resi in range(min_resi, max_resi + 1):
                    f = float(residue_values.get(resi, 0)) / 100.0
                    values.append(f)

                if histone not in results:
                    results[histone] = {}
                if (pdb, chain, offset) not in results[histone]:
                    results[histone][(pdb, chain, offset)] = np.array(values)
                else:
                    # print(results[histone][(pdb, chain, offset)])
                    # print(np.array(values))
                    results[histone][(pdb, chain, offset)] += np.array(values)

    print(name_mapping)

    for histone in results:
        alignment_fname =  d + "interactome/" + ALIGNMENT_FILE.format(histone)
        results_fname =  d + "interactome/" + RESULTS_FILE.format(histone)
        with open(results_fname, 'w') as o, open(alignment_fname, 'w') as fasta:
            o.write("JALVIEW_ANNOTATION\n\n")
            for structure, residues in results[histone].items():
                if residues.sum() == 0.0:
                    continue
                (pdb, chain, offset) = structure
                residues /= residues.sum()
                residues *= 100.0
                residues = residues.round(2)
                values = map("{}".format, residues.tolist())

                if (pdb.upper(), chain) in name_mapping:
                    name = name_mapping[(pdb.upper(), chain)]
                    name2 = name
                else:
                    name = "PDB|{}|{}|{}".format(pdb.upper(), pdb.upper(), chain)
                    name2 = "{}{}".format(pdb.upper(), chain)

                offset = 0
                fasta.write(">{}\n".format(name))
                sequence = sequences[(pdb, chain)]
                fasta.write("\n".join(textwrap.wrap(sequence, 60)) + "\n")

                # Atom(chain=chain, chain_real=chain_real, chain_author=chain_author, resn=resn, resn_short=resn_short, resi=resi, atomn=atomn, atomi=atomi, element=element, xyz=xyz)
                o.write("SEQUENCE_REF\t{}\t{}\n".format(name, offset))
                o.write("BAR_GRAPH\tContacts\tContacts in {}\t{}\n\n".format(name2, "|".join(values) + '|'))
                # BAR_GRAPH 0.0,-|42.857143,K,K 43%|85.71429,R,R 86%; K 14%|


if __name__ == '__main__':
    main()
