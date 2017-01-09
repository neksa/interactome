from flask import Flask, jsonify, abort, Blueprint
from flask import render_template
# from flask import request, redirect
# from flask_restful import Resource, Api

import contextlib
import json

# import urllib
# from random import randint, choice

from urllib.request import urlopen, Request
from urllib.parse import urlencode
from collections import OrderedDict

from sqlalchemy import *
from sqlalchemy.orm import mapper, sessionmaker
from sqlalchemy.ext.automap import automap_base

from collections import defaultdict
import xml.etree.ElementTree as ET

rest_api = Blueprint('rest_api', __name__)

# db_ibis = None
# class ObsInt(object):
#     pass

INT_PPI = 3
INT_SMI = 6
INT_PEPT = 7

REST_PREFIX = "/rest/v1/"

UNIPROT_RESOURCE = "http://www.ebi.ac.uk/proteins/api"
PDB_RESOURCE = "http://www.rcsb.org/pdb/rest"
PUBCHEM_RESOURCE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"

MMDB_RESOURCE = "https://www.ncbi.nlm.nih.gov/Structure/mmdb/mmdb_strview.cgi?uid=27242&format=json"  # &intrac=3


def connect_db_ibis():
    global db_ibis, ObsInt, SidCidInfoHet

    # CONNECT = "mssql+pymssql://anyone:allowed@DDDSQL608/Intrac"
    CONNECT = "mssql+pymssql://anyone:allowed@localhost/Intrac"
    engine = create_engine(CONNECT)
    meta = MetaData(bind=engine)
    meta.reflect(only=['ObsInt', 'PdbXtal', 'MolResFace', 'ExcludedSid', 'IntResFace', 'SidCidInfo', 'StructDomSfam', 'SidCidInfoHet', 'TblMap'])
    Base = automap_base(bind=engine, metadata=meta)
    Base.prepare()

    SidCidInfoHetTable = Table(
        'SidCidInfoHet',
        meta,
        PrimaryKeyConstraint('sid'),
        extend_existing=True)

    ObsInt = Base.classes.ObsInt

    class SidCidInfoHetClass():
        pass

    SidCidInfoHet = mapper(SidCidInfoHetClass, SidCidInfoHetTable)
    # SidCidInfoHet = Base.classes.SidCidInfoHet

    db_ibis = sessionmaker()(bind=engine)


def connect_db_mmdb():
    global db_mmdb, StPdbMap

    # CONNECT = "mssql+pymssql://anyone:allowed@DDDSQL608/PubStructMain"
    CONNECT = "mssql+pymssql://anyone:allowed@localhost/PubStructMain"
    engine = create_engine(CONNECT)
    meta = MetaData(bind=engine)
    meta.reflect(only=['StPdbMap', 'StStructBioUnit', 'StStruct', 'StSid', 'StSeqAccn', 'StBiounitBiopolymers', 'StAsuBiopolymerChain', 'StBiounitLigands'])
    Base = automap_base(bind=engine, metadata=meta)
    Base.prepare()

    StPdbMap = Base.classes.StPdbMap
    db_mmdb = sessionmaker()(bind=engine)


def create_app(conf=None):
    app = Flask(__name__)
    if conf is not None:
        app.config.from_object(conf)
    app.register_blueprint(rest_api)
    connect_db_ibis()
    connect_db_mmdb()
    return app

# app = Flask(__name__)
# api = Api(app)

# CONNECT = "mssql+pymssql://anyone:allowed@MUTAGENE_LOAD/mutagene"


# recs = db.query(ObsInt).filter(ObsInt.mmdb_id == 107035)
# for r in recs:
#     print(r.obs_int_id, r.mol_superfam_acc)


"""
https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/53789435/synonyms/TXT
https://www.ncbi.nlm.nih.gov/Structure/mmdb/mmdb_strview.cgi?uid=27242&format=json
http://www.rcsb.org/pdb/rest/describeMol?structureId=3I5F


curl 'http://www.cancerrxgene.org/translation/drug_list?list=all&export=json&sEcho=1&iColumns=6&sColumns=&iDisplayStart=0&iDisplayLength=30&mDataProp_0=0&mDataProp_1=1&mDataProp_2=2&mDataProp_3=3&mDataProp_4=4&mDataProp_5=5&sSearch=&bRegex=false&sSearch_0=&bRegex_0=false&bSearchable_0=true&sSearch_1=&bRegex_1=false&bSearchable_1=true&sSearch_2=&bRegex_2=false&bSearchable_2=true&sSearch_3=&bRegex_3=false&bSearchable_3=true&sSearch_4=&bRegex_4=false&bSearchable_4=true&sSearch_5=&bRegex_5=false&bSearchable_5=true&iSortCol_0=0&sSortDir_0=asc&iSortingCols=1&bSortable_0=true&bSortable_1=true&bSortable_2=true&bSortable_3=true&bSortable_4=true&bSortable_5=true' \
-XGET \
-H 'Accept: application/json, text/javascript, */*; q=0.01' \
-H 'Referer: http://www.cancerrxgene.org/translation/Drug' \
-H 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10_12_1) AppleWebKit/602.2.14 (KHTML, like Gecko) Version/10.0.1 Safari/602.2.14' \
-H 'X-Requested-With: XMLHttpRequest'


{
    "iTotalRecords": "265",
    "iTotalDisplayRecords": "265",
    "aaData": [
        [
            "<a href=\"http://www.cancerrxgene.org/translation/Drug/1242\">(5Z)-7-Oxozeaenol</a>",
            null,
            "TAK1 (MAP3K7)",
            "other",
            null,
            "951"
        ],




curl 'http://www.cancerrxgene.org/translation/drug_list?list=all&export=json&sEcho=6&iColumns=6&sColumns=&iDisplayStart=0&iDisplayLength=30&mDataProp_0=0&mDataProp_1=1&mDataProp_2=2&mDataProp_3=3&mDataProp_4=4&mDataProp_5=5&sSearch=11433190&bRegex=false&sSearch_0=&bRegex_0=false&bSearchable_0=true&sSearch_1=&bRegex_1=false&bSearchable_1=true&sSearch_2=&bRegex_2=false&bSearchable_2=true&sSearch_3=&bRegex_3=false&bSearchable_3=true&sSearch_4=&bRegex_4=false&bSearchable_4=true&sSearch_5=&bRegex_5=false&bSearchable_5=true&iSortCol_0=0&sSortDir_0=asc&iSortingCols=1&bSortable_0=true&bSortable_1=true&bSortable_2=true&bSortable_3=true&bSortable_4=true&bSortable_5=true' \
-XGET \
-H 'Accept: application/json, text/javascript, */*; q=0.01' \
-H 'Referer: http://www.cancerrxgene.org/translation/Drug' \
-H 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10_12_1) AppleWebKit/602.2.14 (KHTML, like Gecko) Version/10.0.1 Safari/602.2.14' \
-H 'X-Requested-With: XMLHttpRequest'
"""



"""
GETREST_PREFIX +  /structures/protein/<PROTEIN NAME>
Retrieves a list of PDB IDS and chain IDS given a UNIPROT protein name

GETREST_PREFIX +  /structures/compound/<COMPOUND>
Retrieves a list of PDB IDS containing a given PUBCHEM compound

GETREST_PREFIX +  /binding-sites/observed/protein-protein/pdb/<PDB ID>[/<CHAIN>]
Retrieve a list of observed protein-protein binding site residues given a PDB ID (and an optional chain ID).
Each binding site specifies two interacting chains, the interacting domains (when known) and a list of residues with PDB numbering.
Each residue is assigned a conservation score (highly conserved, conserved and not conserved – according to IBIS).

GETREST_PREFIX +  /binding-sites/observed/protein-compound/pdb/<PDB ID>[/<CHAIN>]
Retrieve a list of observed protein-small molecule (PUBCHEM compound) binding site residues given a PDB ID (and an optional chain ID).

GETREST_PREFIX +  /binding-sites/inferred/protein-protein/pdb/<PDB ID>[/<CHAIN>]
Retrieve a list of inferred protein-protein binding site residues given a PDB ID (and an optional chain ID)

-------------------------------

GETREST_PREFIX +  /binding-sites/inferred/protein-compound/pdb/<PDB ID>[/<CHAIN>]
Retrieve a list of inferred protein-small molecule binding site residues given a PDB ID (and an optional chain ID)

GETREST_PREFIX +  /binding-sites/observed/protein-compound/compound/<COMPOUND>
Retrieve a list of observed protein-small molecule binding sites given a small molecule (PUBCHEM compound).

GETREST_PREFIX +  /compounds/similar/compound/<COMPOUND>
Retrieve a list of similar small molecule compounds given a small molecule (PUBCHEM compound ID). Returns MAX 50 items.

GETREST_PREFIX +  /compounds/interacting/protein/<PROTEIN NAME>
Retrieve a list of similar small molecule compounds given a small molecule (PUBCHEM compound ID).

GETREST_PREFIX +  /variants/pdb/<PDB ID>/<CHAIN>/<RESIDUE>
Retrieve all known variants given a residue number in a PDB structure.

"""


def get_mmdb_by_pdb(pdb, chain=None):
    global db_mmdb, StPdbMap

    mmdb_id = 0
    chain_id = 0

    m = db_mmdb.query(StPdbMap).filter(StPdbMap.pdbId == pdb).one()
    mmdb_id = m.mmdbId

    return mmdb_id, chain_id


def get_pdb_by_mmdb(mmdb_id, chain_id=None):
    global db_mmdb, StPdbMap
    pdb = ""
    chain = ""

    m = db_mmdb.query(StPdbMap).filter(StPdbMap.mmdbId == mmdb_id).one()
    pdb = m.pdbId

    # SELECT TOP 1000 [mmdbId]
    #       ,[pdbId]
    #       ,[pdbcls]
    #       ,[expmthd]
    #       ,[reso]
    #       ,[ecno]
    #       ,[rlsdate]
    #       ,[depdate]
    #   FROM [PubStructMain].[dbo].[StPdbMap]
    return pdb, chain


# def get_observed_bs(int_type, mmdb_id, chain_id=None):
#     pass


# def get_inferred_bs(int_type, mmdb_id, chain_id=None):
#     pass


@rest_api.route(REST_PREFIX + "structures/protein/<protein_name>")
def get_structures_by_protein_name(protein_name):
    """
    Uniprot API - get PDB list
    """
    # pdb&format=tab&compress=yes&columns=id,database(PDB)
    # urlparams = urlencode({'query': '"' + protein_name + '" AND organism:9606 AND database:pdb', 'format': 'tab', 'columns': 'id,database(PDB)'})
    urlparams = urlencode({'gene': protein_name, 'taxid': '9606', 'format': 'tab', 'columns': 'id'})
    url = UNIPROT_RESOURCE + "/proteins?%s" % urlparams
    print(url)
    r = Request(url, headers={'Accept': 'text/x-fasta'})
    accessions = set()
    try:
        # http://www.ncbi.nlm.nih.gov/Structure/ibis/ibis.cgi?cmd=text&id=734936&type=3&ofm=3
        with contextlib.closing(urlopen(r)) as f:
            tab = f.read().decode('utf-8')
            for line in tab.split("\n"):
                if line.startswith(">"):
                    cols = line[1:].split("|")
                    if cols[0] != 'sp':
                        continue
                    accession = cols[1]
                    accessions.add(accession)
    except:
        return abort(404)

    print(accessions)

    pdb_query = """
    <orgPdbQuery>
        <queryType>org.pdb.query.simple.UpAccessionIdQuery</queryType>
        <description></description>
        <accessionIdList>{}</accessionIdList>
    </orgPdbQuery>
    """.format(",".join(accessions))

    print(pdb_query)
    urlparams = urlencode({'sortfield': 'Resolution'})
    url = PDB_RESOURCE + "/search?%s" % urlparams
    print(url)
    pdbs = set()
    r = Request(url, data=pdb_query.encode('ascii'), method='POST')
    try:
        with contextlib.closing(urlopen(r)) as f:
            tab = f.read().decode('utf-8')
            for line in tab.split("\n"):
                if ":" not in line:
                    continue
                pdb, entity = line.split(":")
                pdbs.add((pdb, int(entity)))
    except:
        raise
        return abort(404)

    pdb_chains = defaultdict(list)

    query_ids = ",".join(list(set([pdb[0].upper() for pdb in pdbs])))
    URI = "http://www.rcsb.org/pdb/rest/describeMol?structureId=" + query_ids
    r = Request(URI)
    try:
        with urlopen(r) as f:
            text = f.read().decode('utf-8')
            root = ET.fromstring(text)
            structures = root.findall("./structureId")
            for struct in structures:
                pdb = struct.attrib['id']
                polymers = struct.findall('./polymer')
                for polymer in polymers:
                    entity = polymer.attrib['entityNr']
                    chains = polymer.findall('./chain')
                    for ch in chains:
                        c = ch.attrib['id']
                        pdb_chains[(pdb, int(entity))].append(c)
    except:
        raise

    response = []
    for pdb, entity in pdbs:
        chains = pdb_chains[(pdb, entity)]
        if len(chains) > 0:
            for c in chains:
                response.append({
                    'pdb': pdb,
                    'chain': c})

    return jsonify({'response': response})


@rest_api.route(REST_PREFIX + "structures/compound/<compound>")
def get_structures_by_compound(compound):
    """
    Pubchem API
    """

    url = PUBCHEM_RESOURCE + "compound/cid/{}/xrefs/MMDBID/JSON".format(compound)  # 2244
    try:
        with contextlib.closing(urlopen(url)) as f:
            result = f.read().decode('utf-8')
            info = json.loads(result)
            mmdbs = []
            for rec in info['InformationList']['Information']:
                for mmdb in rec['MMDBID']:
                    mmdbs.append(mmdb)
            pdbs = []
            for m in mmdbs:
                pdb, chain = get_pdb_by_mmdb(int(m))
                pdbs.append(pdb)

            return jsonify({'response': pdbs})

    except:
        abort(404)
    return jsonify({'response': ""})


@rest_api.route(REST_PREFIX + "compounds/similar/compound/<compound>")
def get_compounds_by_similar_compound(compound):
    """
    Pubchem API
    """

    # Similar 3D
    cids_3d = []

    # url = PUBCHEM_RESOURCE + "compound/fastsimilarity_3d/cid/{}/cids/JSON".format(compound)  # 2244
    # print(url)
    # try:
    #     with contextlib.closing(urlopen(url)) as f:
    #         result = f.read().decode('utf-8')
    #         info = json.loads(result)
    #         for cid in info['IdentifierList']["CID"]:
    #             cids_3d.append(cid)
    # except:
    #     raise
    #     abort(500)

    # print(cids_3d)

    # Similar 2D
    url = PUBCHEM_RESOURCE + "compound/fastsimilarity_2d/cid/{}/cids/JSON".format(compound)  # 2244
    print(url)
    cids_2d = []
    try:
        with contextlib.closing(urlopen(url)) as f:
            result = f.read().decode('utf-8')
            info = json.loads(result)
            for cid in info['IdentifierList']["CID"]:
                cids_2d.append(cid)
    except:
        raise
        abort(500)

    print(cids_2d)

    # cids = list(set(cids_3d) & set(cids_2d))
    cids = list(set(cids_3d) | set(cids_2d))
    if len(cids) > 0:
        return jsonify({'response': cids})
    else:
        return jsonify({'error': "Not found"})


@rest_api.route(REST_PREFIX + "binding-sites/<obs_inf>/<int_type>/pdb/<pdb>")
@rest_api.route(REST_PREFIX + "binding-sites/<obs_inf>/<int_type>/pdb/<pdb>/<chain>")
def get_bs_by_pdb(obs_inf, int_type, pdb, chain=None):
    """
    PDB to MMDB
    IBIS query observed PPI
    """
    print(obs_inf)
    print(int_type)

    if obs_inf not in ("observed", "inferred"):
        print("11")
        abort(404)

    if int_type not in ("protein-protein", "protein-compound", "protein-peptide"):
        print("22222")
        abort(404)

    inferred = obs_inf == "inferred"

    interaction = INT_PPI if int_type == "protein-protein" else False
    interaction = INT_SMI if int_type == "protein-compound" else False
    interaction = INT_PEPT if int_type == "protein-peptide" else False

    get_obs_smi_int(pdb)
    return ibis(inferred, interaction, pdb, chain)
    # return jsonify({'response': ""})


def get_obs_smi_int(pdb):
    global db_ibis, ObsInt, SidCidInfoHet

    mmdb_id, chain = get_mmdb_by_pdb(pdb)
    q = db_ibis.query(ObsInt).filter(ObsInt.mmdb_id == mmdb_id).all()
    data = []
    for a in q:
        sid = a.int_sid
        cid = db_ibis.query(SidCidInfoHet).get(sid)
        print(sid, cid)
        data.append(a.__dict__)

    import pprint
    pprint.pprint(data)


def test():
    global db_ibis, ObsInt, SidCidInfoHet

    mmdb_id, chain = get_mmdb_by_pdb(pdb)
    q = db_ibid.query(ObsInt).filter(ObsInt.mmdb_id == mmdb_id).all()
    data = []
    for a in q:
        sid = a.int_sid
        cid = db_ibid.query(SidCidInfoHet).get(sid)
        print(sid, cid)
        data.append(a.__dict__)

    import pprint
    pprint.pprint(data)

    # [obs_int_id]
    #      ,[mmdb_id]
    #      ,[mol_sdi_id]
    #      ,[mol_gi]
    #      ,[mol_mol_id]
    #      ,[type]
    #      ,[image_id]
    #      ,[mol_superfam_id]
    #      ,[mol_superfam_acc]
    #      ,[mol_cd_from]
    #      ,[mol_cd_to]
    #      ,[mol_cd_pssmid]
    #      ,[mol_taxid]
    #      ,[int_mol_id]
    #      ,[int_gi]
    #      ,[int_sid]
    #      ,[mol_cd_cur_annot]
    #      ,[int_sdi_id]
    #      ,[has_int_resface]
    #      ,[int_superfam_id]
    #      ,[int_superfam_acc]
    #      ,[int_cd_from]
    #      ,[int_cd_to]
    #      ,[int_cd_pssmid]
    #      ,[int_taxid]
    #      ,[n_resface_contacts]
    #      ,[mol_cd_cur_pssm]
    #      ,[int_sequence]
    #      ,[molface_counts]
    #      ,[pisa_status]




# @rest_api.route(REST_PREFIX + "binding-sites/observed/protein-compound/pdb/<pdb>")
# @rest_api.route(REST_PREFIX + "binding-sites/observed/protein-compound/pdb/<pdb>/<chain>")
# def get_observed_smi_bs_by_pdb(pdb, chain=None):
#     """
#     PDB to MMDB
#     IBIS query observed SMI
#     """
#     # return ibis(pdb, chain)
#     return jsonify({'response': ""})


# @rest_api.route(REST_PREFIX + "binding-sites/inferred/protein-protein/pdb/<pdb>")
# @rest_api.route(REST_PREFIX + "binding-sites/inferred/protein-protein/pdb/<pdb>/<chain>")
# def get_inferred_ppi_bs_by_pdb(pdb, chain=None):
#     """
#     PDB to MMDB
#     IBIS query observed PPI
#     """
#     return ibis(pdb, chain)
#     return jsonify({'response': ""})


# @rest_api.route(REST_PREFIX + "binding-sites/inferred/protein-compound/pdb/<pdb>")
# @rest_api.route(REST_PREFIX + "binding-sites/inferred/protein-compound/pdb/<pdb>/<chain>")
# def get_inferred_smi_bs_by_pdb(pdb, chain=None):
#     """
#     PDB to MMDB
#     IBIS query observed SMI
#     """
#     # return ibis(pdb, chain)
#     return jsonify({'response': ""})


# @rest_api.route(REST_PREFIX + "binding-sites/observed/protein-compound/compound/<compound>")
# def get_observed_bs_by_compound(compound):
#     """
#     IBIS query observed SMI by compound
#     """
#     return jsonify({'response': ""})


@rest_api.route(REST_PREFIX + "compounds/interacting/protein/<protein_name>")
def get_compounds_by_protein(protein_name):
    """
    UNIPROT API query pdb by protein_name
    IBIS convert pdb to MMDB
    IBIS query observed SMI by MMDB
    """
    return jsonify({'response': ""})


@rest_api.route(REST_PREFIX + "variants/pdb/<pdb>/<chain>/<residue>")
def get_variants_by_pdb(pdb, chain, residue):
    """
    dbSNP query variants by PDB ?????
    """
    return jsonify({'response': ""})


@rest_api.route("/")
def home():
    # recs = db_ibis.query(ObsInt).filter(ObsInt.mmdb_id == 107035)
    # for r in recs:
    #     print(r.obs_int_id, r.mol_superfam_acc)

    return render_template('rest.html', prefix=REST_PREFIX)
    # return("REST API for drug discovery pipeline based on IBIS and PubChem")


# @app.route("/ibis/<pdb>/<chain>")
def ibis(inferred, interaction, pdb, chain):
    chain = chain if chain else ""

    sdi = None
    other_sdis = []
    interac = []

    # standard_error = " "
    params = {'pdb': pdb, 'chain': chain}

    # Retrieve query SDI:
    try:
        urlparams = urlencode({'cmd': 'sumtbl', 'search': pdb + chain})
        url = "https://www.ncbi.nlm.nih.gov/Structure/ibis/ibis.cgi?%s" % urlparams
        with contextlib.closing(urlopen(url)) as f:
            html = f.read().decode('utf-8')
            # return html
            for line in html.split("\n"):
                if "Protein-protein" in line:
                    for k in line.split():
                        if k.startswith("sdi"):
                            sdi = k[5:-1]
                            sdi = int(sdi)
        if sdi is None:
            print("No SDI")
            abort(404)
    except:
        raise
        abort(404)

    try:
        urlparams = urlencode({'cmd': 'text', 'type': 3, 'id': sdi, 'ofm': 3})
        url = "https://www.ncbi.nlm.nih.gov/Structure/ibis/ibis.cgi?%s" % urlparams
        # print(url)
        # http://www.ncbi.nlm.nih.gov/Structure/ibis/ibis.cgi?cmd=text&id=734936&type=3&ofm=3
        with contextlib.closing(urlopen(url)) as f:
            csv = f.read().decode('utf-8')

            for j, line in enumerate(csv.split("\n")):
                if j == 0 or line == "":
                    continue
                if line.startswith("<html"):
                    break
                fields = line.split(",")
                # fields[3]  # MMDB
                resi = fields[4].strip().split()  # PDB
                resn = fields[5].strip().split()  # RES
                conserv = fields[6].strip().split()  # conservation

                singleton = fields[7].strip()  # singleton
                domain = fields[-2]  # DOMAIN

                site = []
                for i, x in enumerate(resi):
                    site.append({
                        'i': i,
                        'resi': resi[i],
                        'resn': resn[i],
                        'conserv': conserv[i],
                        })
                interac.append({
                    'partner': domain,
                    'site': site,
                    'singleton': singleton,
                    })

        if len(interac) == 0:
            print("EMPTY INTERAC")
            abort(404)
    except:
        raise
        abort(404)
    try:
        urlparams = urlencode({'type': 3, 'id': sdi})
        url = "https://www.ncbi.nlm.nih.gov/Structure/ibis/ibis.cgi?%s" % urlparams
        with contextlib.closing(urlopen(url)) as f:
            html = f.read().decode('utf-8')
            for line in html.split("\n"):
                line = line.strip()
                if line.startswith("InfId2Info"):
                    json_string = line.split("=", 1)[-1][:-1].replace("'", "\"")
                    # print("JSON: ", json_string)
                    clusters = json.loads(json_string, object_pairs_hook=OrderedDict)

                if line.startswith("var FaceResNumb"):
                    json_string = line.split("=", 1)[-1][:-1].replace("'", "\"").replace("],", "],\"").replace("{", "{\"").replace(":[", "\":[")
                    # print("JSON: ", json_string)
                    cluster_resi = json.loads(json_string, object_pairs_hook=OrderedDict)

                # if line.startswith("var FaceResColor"):
                #     json_string = line.split("=", 1)[-1][:-1].replace("'", "\"").replace("],", "],\"").replace("{", "{\"").replace(":[", "\":[")
                #     # print("JSON: ", json_string)
                #     cluster_color = json.loads(json_string, object_pairs_hook=OrderedDict)

                # AJAX: <div align="center" class="seqgph" id="seqgph"><img alt="Aligned Interaction Positions" gid="&amp;id=326283&amp;type=3&amp;focus=452075" id="bigloader" src="images/bigLoader.gif"></div>
                # if line.startswith("<area coords") and "conserved domain information" in line:
                #     print(line)
                #     for subitem in line.split():
                #         if "uid" in subitem:
                #             uid = int(subitem.split("=")[-1][:-1])
                #             print("Identified UID", uid)
                #             if uid != sdi:
                #                 other_sdis.append(uid)

        cluster_codes = list(clusters.keys())
        for j, x in enumerate(interac):
            interac[j]['homologs'] = []
            code = cluster_codes[j]
            # print("HOMOLOGS:")
            for member in clusters[code][1:-1]:
                member_pdb = member[0]
                member_chain1, _, member_chain2, _ = member[1].split("::", 3)
                member_seq = member[4]

                site = []
                for i, s in enumerate(interac[j]['site']):
                    site.append({
                        'resi': int(cluster_resi[code][i]),
                        'resn': member_seq[i],
                        'conserv': int(s['conserv']),
                        })
                interac[j]['homologs'].append({
                    'pdb': member_pdb,
                    'chain1': member_chain1,
                    'chain2': member_chain2,
                    'site': site
                })
    except:
        raise
        abort(404)

    params['interac'] = interac
    params['other_sdis'] = other_sdis
    params['sdi'] = sdi
    return jsonify({'response': params})


if __name__ == "__main__":
    app = create_app()
    app.run(debug=True)


    """
    /****** Script for SelectTopNRows command from SSMS  ******/
    SELECT TOP 1000 [obs_int_id]
          ,[mmdb_id]
          ,[mol_sdi_id]
          ,[mol_gi]
          ,[mol_mol_id]
          ,[type]
          ,[image_id]
          ,[mol_superfam_id]
          ,[mol_superfam_acc]
          ,[mol_cd_from]
          ,[mol_cd_to]
          ,[mol_cd_pssmid]
          ,[mol_taxid]
          ,[int_mol_id]
          ,[int_gi]
          ,[int_sid]
          ,[mol_cd_cur_annot]
          ,[int_sdi_id]
          ,[has_int_resface]
          ,[int_superfam_id]
          ,[int_superfam_acc]
          ,[int_cd_from]
          ,[int_cd_to]
          ,[int_cd_pssmid]
          ,[int_taxid]
          ,[n_resface_contacts]
          ,[mol_cd_cur_pssm]
          ,[int_sequence]
          ,[molface_counts]
          ,[pisa_status]
      FROM [Intrac].[dbo].[ObsInt]



    /****** Script for SelectTopNRows command from SSMS  ******/
    SELECT TOP 1000 [sid]
          ,[cid]
          ,[sidname]
          ,[cidname]
          ,[exclude]
          ,[mesh]
          ,[pharmaction]
          ,[active_assays]
          ,[assays]
          ,[pdbhet]
          ,[mmcif]
      FROM [Intrac].[dbo].[SidCidInfoHet]
    """



"""
/****** Script for SelectTopNRows command from SSMS  ******/
SELECT TOP 1000 [pdbId]
      ,[chnLett]
      ,[acxn]
      ,[gi]
      ,[mmdbId]
      ,[geneId]
      ,[molId]
      ,[pig]
  FROM [PubStructMain].[dbo].[StSeqAccn]


  /****** Script for SelectTopNRows command from SSMS  ******/
SELECT TOP 1000 [mmdbId]
      ,[molId]
      ,[lignme]
      ,[sid]
      ,[liglngnme]
      ,[ligsyn]
  FROM [PubStructMain].[dbo].[StSid]


/****** Script for SelectTopNRows command from SSMS  ******/
SELECT TOP 1000 [mmdbId]
      ,[pdbId]
      ,[pdbcls]
      ,[expmthd]
      ,[reso]
      ,[ecno]
      ,[rlsdate]
      ,[depdate]
  FROM [PubStructMain].[dbo].[StPdbMap]


/****** Script for SelectTopNRows command from SSMS  ******/
SELECT TOP 1000 [acc]
      ,[mmdbId]
      ,[asuMolId]
      ,[pdbMolId]
      ,[chnLett]
      ,[kind]
      ,[chnLettPrefix]
  FROM [PubStructMain].[dbo].[StAsuBiopolymerChain]

/****** Script for SelectTopNRows command from SSMS  ******/
SELECT TOP 1000 [buAcc]
      ,[molId]
      ,[asuChainAcc]
      ,[chainNme]
  FROM [PubStructMain].[dbo].[StBiounitBiopolymers]

  /****** Script for SelectTopNRows command from SSMS  ******/
SELECT TOP 1000 [buAcc]
      ,[molId]
      ,[chemMolId]
      ,[kind]
  FROM [PubStructMain].[dbo].[StBiounitLigands]


"""



"""
<molDescription>
<structureId id="4HHB">
<polymer entityNr="1" length="141" type="protein" weight="15150.5">
<chain id="A"/>
<chain id="C"/>
<Taxonomy name="Homo sapiens" id="9606"/>
<macroMolecule name="Hemoglobin subunit alpha">
<accession id="P69905"/>
</macroMolecule>
<polymerDescription description="HEMOGLOBIN (DEOXY) (ALPHA CHAIN)"/>
</polymer>
<polymer entityNr="2" length="146" type="protein" weight="15890.4">
<chain id="B"/>
<chain id="D"/>
<Taxonomy name="Homo sapiens" id="9606"/>
<macroMolecule name="Hemoglobin subunit beta">
<accession id="P68871"/>
</macroMolecule>
<polymerDescription description="HEMOGLOBIN (DEOXY) (BETA CHAIN)"/>
</polymer>
</structureId>
</molDescription>
"""


