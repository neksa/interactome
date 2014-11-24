"""
Load interactions to database
"""

from sqlalchemy import create_engine
from sqlalchemy import Table, Column, Integer, Float, String, MetaData, ForeignKey
# from sqlalchemy.sql import text
from itertools import islice
engine = create_engine('postgresql://localhost/agoncear')#, echo=True)

metadata = MetaData()
species = Table('species', metadata,
                Column('tax_id', Integer, primary_key=True),
                Column('name', String))

pairs =   Table('pairs', metadata,
                Column('pair_id', Integer, primary_key=True),
                Column('tax_id', Integer, ForeignKey('species.tax_id'), nullable=False),
                Column('rev', Integer, nullable=False),
                Column('id_a', String, nullable=False),
                Column('id_b', String, nullable=False),
                Column('uniprot_ac_a', String, nullable=False),
                Column('uniprot_ac_b', String, nullable=False),
                Column('tpl', String),
                Column('query_type', String),
                Column('template_type', String),
                Column('score_model_z', Float),
                Column('score_model_minus_avg', Float),
                Column('score3', Float),
                Column('score4', Float),
                Column('score5', Float),
                Column('score6', Float),
                Column('identical_a', Integer),
                Column('identical_b', Integer),
                Column('positive_a', Integer),
                Column('positive_b', Integer),
                Column('aln_len_a', Integer),
                Column('aln_len_b', Integer),
                Column('bs_len_a', Integer),
                Column('bs_len_b', Integer),
                Column('bs_covered_a', Integer),
                Column('bs_covered_b', Integer),
                Column('bs_aligned_a', Integer),
                Column('bs_aligned_b', Integer),
                Column('bs_identical_a', Integer),
                Column('bs_identical_b', Integer),
                Column('bs_positive_a', Integer),
                Column('bs_positive_b', Integer),
                Column('bs_contacts_a', Integer),
                Column('bs_contacts_b', Integer),
                Column('bs_blosum_a', Float),
                Column('bs_blosum_b', Float),
                Column('bs_score1_a', Float),
                Column('bs_score1_b', Float),
                Column('site_a', String),
                Column('site_b', String))

proteins = Table('proteins', metadata,
                 Column('uniprot_ac', String, primary_key=True),
                 Column('uniprot_accessions', String),
                 Column('uniprot_id', String),
                 Column('tax_id', Integer, ForeignKey('species.tax_id'), nullable=False),
                 Column('name', String),
                 Column('sequence', String),
                 Column('seqlen', Integer))

idmapping = Table('idmapping', metadata,
                  Column('uniprot_ac', String, primary_key=True, nullable=True),
                  Column('uniprot_id', String),
                  Column('gene_id', Integer),
                  Column('gi', Integer),
                  Column('refseq', String),
                  Column('pdb', String),
                  Column('go', String),
                  Column('ensembl', String))

metadata.create_all(engine)
con = engine.connect()


def setup_database():
    pass


def load_uniprot():
    pass


def load_interactions():

    def load_gi_uniprot_mapping(fname):
        mapping = {}
        inv_mapping = {}
        with open(fname) as f:
            for line in f:
                gis, uniprot = line.strip().split()
                for gi in gis.split(","):
                    mapping[gi] = uniprot
                    inv_mapping[uniprot] = gi
        return mapping, inv_mapping

    data_dir = "/Users/agoncear/data/Ecoli/"
    mapping_fname = data_dir + "gi_uniprot_mapping.tab"
    matches_fname = "/Users/agoncear/projects/Interactome/Workflow/Alignments/matches_ecoli.tab"
    gi_mapping, uniprot_mapping = load_gi_uniprot_mapping(mapping_fname)

    with open(matches_fname) as f:
        for line in islice(f, 1, None):
            # print "L", len(line.strip().split("\t"))
            # for i, a in enumerate(line.strip().split("\t")):
            #     print i, "==", a, "=="
            # 11 + 2 + 11 + 11 = 35
            try:
                queryA, queryB, tpl, query_type, template_type, \
                    SC1, SC2, SC3, SC4, SC5, SC6, \
                    identicalA, positiveA, aln_lenA, bs_lenA, bs_coveredA, bs_alignedA, bs_identicalA, bs_positiveA, bs_contactsA, bs_BLOSUMA, bs_score1A, \
                    identicalB, positiveB, aln_lenB, bs_lenB, bs_coveredB, bs_alignedB, bs_identicalB, bs_positiveB, bs_contactsB, bs_BLOSUMB, bs_score1B, \
                    siteA, siteB = line.strip().split("\t")
                    # _, _, _, _, \

                score_template_full, score_template, score_model, scaled_score_template_full, scaled_score_template, scaled_score = map(
                    float, (SC1, SC2, SC3, SC4, SC5, SC6))

                zscore, model_minus_avg = score_template_full, score_template

                identicalA, positiveA, aln_lenA, bs_lenA, bs_coveredA, bs_alignedA, bs_identicalA, bs_positiveA, bs_contactsA, bs_BLOSUMA, bs_score1A = map(
                    float, (identicalA, positiveA, aln_lenA, bs_lenA, bs_coveredA, bs_alignedA, bs_identicalA, bs_positiveA, bs_contactsA, bs_BLOSUMA, bs_score1A))

                identicalB, positiveB, aln_lenB, bs_lenB, bs_coveredB, bs_alignedB, bs_identicalB, bs_positiveB, bs_contactsB, bs_BLOSUMB, bs_score1B = map(
                    float, (identicalB, positiveB, aln_lenB, bs_lenB, bs_coveredB, bs_alignedB, bs_identicalB, bs_positiveB, bs_contactsB, bs_BLOSUMB, bs_score1B))

                # YMR101C SRT1 SGDID:S000004707, Chr XIII from....
                queryA = queryA.split()[0]
                queryB = queryB.split()[0]

                uniprotA = gi_mapping.get(queryA)
                uniprotB = gi_mapping.get(queryB)

                if uniprotA is None or uniprotB is None:
                    continue

                ins = pairs.insert().values(
                    tax_id=511145,
                    rev=1,
                    id_a=queryA,
                    id_b=queryB,
                    uniprot_ac_a=uniprotA,
                    uniprot_ac_b=uniprotB,
                    tpl=tpl,
                    query_type=query_type,
                    template_type=template_type,
                    score_model_z=zscore,
                    score_model_minus_avg=model_minus_avg,
                    score3=score_model,
                    score4=SC4,
                    score5=SC5,
                    score6=SC6,
                    identical_a=identicalA,
                    identical_b=identicalB,
                    positive_a=positiveA,
                    positive_b=positiveB,
                    aln_len_a=aln_lenA,
                    aln_len_b=aln_lenB,
                    bs_len_a=bs_lenA,
                    bs_len_b=bs_lenB,
                    bs_covered_a=bs_coveredA,
                    bs_covered_b=bs_coveredB,
                    bs_aligned_a=bs_alignedA,
                    bs_aligned_b=bs_alignedB,
                    bs_identical_a=bs_identicalA,
                    bs_identical_b=bs_identicalB,
                    bs_positive_a=bs_positiveA,
                    bs_positive_b=bs_positiveB,
                    bs_contacts_a=bs_contactsA,
                    bs_contacts_b=bs_contactsB,
                    bs_blosum_a=bs_BLOSUMA,
                    bs_blosum_b=bs_BLOSUMB,
                    bs_score1_a=bs_score1A,
                    bs_score1_b=bs_score1B,
                    site_a=siteA,
                    site_b=siteB)
                con.execute(ins)

            except:
                print "Error in line", line
                raise

    pass


def load_idmapping():
    pass


def load_species():
    ins = species.insert().values(tax_id=511145, name='Escherichia coli str. K-12 substr. MG1655')
    con.execute(ins)


def main():
    # load_species()
    load_interactions()


if __name__ == '__main__':
    main()
