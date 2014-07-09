class SIFTS:
    def getPDBMapping(self, fname):
        """
        Mapping of PDB:CHAIN to Swissprot, pdb_chain_uniprot.tsv

        # 2014/04/28 - 12:58
        PDB CHAIN   SP_PRIMARY  RES_BEG RES_END PDB_BEG PDB_END SP_BEG  SP_END
        101m    A   P02185  1   154 0   153 1   154
        102l    A   P00720  1   40  1   40  1   40

        1a25    A       P68403  14      149     154     289     154     289
        1a25    B       P68403  14      149     154     289     154     289
        """
        PDBChain = namedtuple('PDBChain', 'protein, seqres_from, seqres_to, atom_from, atom_to, protein_from, protein_to')
        mapping = dict()
        # letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        with open(fname) as f:
            for line in islice(f, 2, None):
                try:
                    pdb, chain, sp, res_from, res_to, pdb_from, pdb_to, sp_from, sp_to = line.strip().split("\t", 8)
                    pdb_from = re.sub(r'[A-Z]', '', pdb_from)
                    pdb_to = re.sub(r'[A-Z]', '', pdb_to)
                    # pdb_from = pdb_from.translate(None, 'ABCDEFGHIJKLMNOPQRSTUVWXYZ')
                    # pdb_to = pdb_to.translate(None, )
                    # mapping[pdb.upper() + chain.upper()] = {
                    mapping[pdb + chain] = PDBChain(protein = sp,
                            seqres_from = int(res_from), seqres_to = int(res_to),
                            atom_from = int(pdb_from), atom_to = int(pdb_to),
                            protein_from = int(sp_from), protein_to = int(sp_to))
                except:
                    print "Error: skipping line in SIFTS mapping", line.strip()
        return mapping
