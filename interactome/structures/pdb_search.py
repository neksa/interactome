

class PDBSearchResults:
    def getTabularResults(self, fname):
        records = []
        PDBSearchRecord = namedtuple('PDBSearchRecord', 'pdb, chain, resolution, source, tax, entity, sequence, chain_len, db_id')
        with open(fname) as f:
            reader = csv.reader(islice(f, 1, None), delimiter=',', quotechar='"')
            for r in reader:
                if len(r) > 0:
                    rec = PDBSearchRecord._make(r)
                    # yield rec
                    records.append(rec)
            return records    
            # for rec in map(PDBSearchRecord._make, reader):
            #     print rec.pdb, rec.chain, rec.resolution
            #     print rec.sequence
            #     print
            #     records.append(rec)

    def filterPDB(self, records):
        excluded_pdb = []
        min_chain_len = 50
        tax_id = 9606
        for rec in records:
            try:
                if int(rec.tax) != tax_id:
                    excluded_pdb.append(rec.pdb)
                    continue

                if int(rec.chain_len) < min_chain_len:
                    excluded_pdb.append(rec.pdb)
                    continue
            except:
                excluded_pdb.append(rec.pdb)
                continue                

        excluded_pdb = set(excluded_pdb)
        new_records = [rec for rec in records if rec.pdb not in excluded_pdb]
        return new_records


    def searchToFasta(self, records, fasta_file):
        wrapper = TextWrapper(break_on_hyphens = False, width = 60)
        with open(fasta_file, 'w') as o:
            for r in records:
                o.write(">{}{}\n".format(r.pdb.upper(), r.chain.upper()))
                o.write("{}\n".format(wrapper.fill(r.sequence)))
