from flask import Flask, jsonify, abort, Blueprint
# from flask_restful import Resource, Api
from flask import request, render_template, redirect

import contextlib
import json
# import urllib

# from random import choice
from random import randint

from urllib.request import urlopen, Request
from urllib.parse import urlencode
from collections import OrderedDict

from sqlalchemy import *
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.automap import automap_base

rest_api = Blueprint('rest_api', __name__)


def connect_db_ibis():
    global ibis_db, ObsInt

    CONNECT = "mssql+pymssql://anyone:allowed@DDDSQL608/Intrac"
    engine = create_engine(CONNECT)
    meta = MetaData(bind=engine)
    meta.reflect(only=['ObsInt', 'PdbXtal', 'MolResFace', 'ExcludedSid', 'IntResFace', 'SidCidInfo', 'StructDomSfam', 'SidCidInfoHet', 'TblMap'])
    Base = automap_base(bind=engine, metadata=meta)
    Base.prepare()

    ObsInt = Base.classes.ObsInt
    ibis_db = sessionmaker()(bind=engine)


def create_app(conf=None):
    app = Flask(__name__)
    if conf is not None:
        app.config.from_object(conf)
    app.register_blueprint(rest_api)
    connect_db_ibis()
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
    mmdb_id = 0
    chain_id = 0
    return mmdb_id, chain_id


def get_pdb_by_mmdb(mmdb_id, chain_id=None):
    pdb = ""
    chain = ""
    return pdb, chain


INT_PPI = 3
INT_SMI = 6

UNIPROT_RESOURCE = "http://www.ebi.ac.uk/proteins/api"
PDB_RESOURCE = "http://www.rcsb.org/pdb/rest"
PUBCHEM_RESOURCE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
REST_PREFIX = "/rest/pipeline"

MMDB_RESOURCE = "https://www.ncbi.nlm.nih.gov/Structure/mmdb/mmdb_strview.cgi?uid=27242&format=json"  # &intrac=3


def get_observed_bs(int_type, mmdb_id, chain_id=None):
    pass


def get_inferred_bs(int_type, mmdb_id, chain_id=None):
    pass


@rest_api.route(REST_PREFIX + "/structures/protein/<protein_name>")
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

    response = []
    for pdb, entity in pdbs:
        response.append({
            'pdb': pdb,
            'entity': int(entity)})

    return jsonify({'response': response})


@rest_api.route(REST_PREFIX + "/structures/compound/<compound>")
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
            return jsonify({'response': mmdbs})
    except:
        abort(404)
    return jsonify({'response': ""})


@rest_api.route(REST_PREFIX + "/compounds/similar/compound/<compound>")
def get_compounds_by_similar_compound(compound):
    """
    Pubchem API
    """

    # url = PUBCHEM_RESOURCE + "compound/fastsimilarity_3d/cid/{}/cids/JSON".format(compound)  # 2244
    url = PUBCHEM_RESOURCE + "compound/fastsimilarity_2d/cid/{}/cids/JSON".format(compound)  # 2244
    # print(url)
    try:
        with contextlib.closing(urlopen(url)) as f:
            result = f.read().decode('utf-8')
            info = json.loads(result)
            cids = []
            for cid in info['IdentifierList']["CID"]:
                cids.append(cid)
            return jsonify({'response': cids})
    except:
        # raise
        abort(404)
    return jsonify({'response': ""})


@rest_api.route(REST_PREFIX + "/binding-sites/observed/protein-protein/pdb/<pdb>")
@rest_api.route(REST_PREFIX + "/binding-sites/observed/protein-protein/pdb/<pdb>/<chain>")
def get_observed_ppi_bs_by_pdb(pdb, chain=None):
    """
    PDB to MMDB
    IBIS query observed PPI
    """
    return ibis(pdb, chain)
    # return jsonify({'response': ""})


@rest_api.route(REST_PREFIX + "/binding-sites/observed/protein-compound/pdb/<pdb>")
@rest_api.route(REST_PREFIX + "/binding-sites/observed/protein-compound/pdb/<pdb>/<chain>")
def get_observed_smi_bs_by_pdb(pdb, chain=None):
    """
    PDB to MMDB
    IBIS query observed SMI
    """
    # return ibis(pdb, chain)
    return jsonify({'response': ""})


@rest_api.route(REST_PREFIX + "/binding-sites/inferred/protein-protein/pdb/<pdb>")
@rest_api.route(REST_PREFIX + "/binding-sites/inferred/protein-protein/pdb/<pdb>/<chain>")
def get_inferred_ppi_bs_by_pdb(pdb, chain=None):
    """
    PDB to MMDB
    IBIS query observed PPI
    """
    return ibis(pdb, chain)
    return jsonify({'response': ""})


@rest_api.route(REST_PREFIX + "/binding-sites/inferred/protein-compound/pdb/<pdb>")
@rest_api.route(REST_PREFIX + "/binding-sites/inferred/protein-compound/pdb/<pdb>/<chain>")
def get_inferred_smi_bs_by_pdb(pdb, chain=None):
    """
    PDB to MMDB
    IBIS query observed SMI
    """
    # return ibis(pdb, chain)
    return jsonify({'response': ""})


@rest_api.route(REST_PREFIX + "/binding-sites/observed/protein-compound/compound/<compound>")
def get_observed_bs_by_compound(compound):
    """
    IBIS query observed SMI by compound
    """
    return jsonify({'response': ""})


@rest_api.route(REST_PREFIX + "/compounds/interacting/protein/<protein_name>")
def get_compounds_by_protein(protein_name):
    """
    UNIPROT API query pdb by protein_name
    IBIS convert pdb to MMDB
    IBIS query observed SMI by MMDB
    """
    return jsonify({'response': ""})


@rest_api.route(REST_PREFIX + "/variants/pdb/<pdb>/<chain>/<residue>")
def get_variants_by_pdb(pdb, chain, residue):
    """
    dbSNP query variants by PDB ?????
    """
    return jsonify({'response': ""})


@rest_api.route("/")
def home():
    return render_template('rest.html')
    # return("REST API for drug discovery pipeline based on IBIS and PubChem")


# @app.route("/ibis/<pdb>/<chain>")
def ibis(pdb, chain):

    sdi = None
    other_sdis = []
    interac = []

    # standard_error = " "
    params = {'pdb': pdb, 'chain': chain}

    try:
        urlparams = urlencode({'cmd': 'sumtbl', 'search': pdb + chain})
        url = "http://www.ncbi.nlm.nih.gov/Structure/ibis/ibis.cgi?%s" % urlparams
        with contextlib.closing(urlopen(url)) as f:
            html = f.read().decode('utf-8')
            for line in html.split("\n"):
                if "Protein-protein" in line:
                    for k in line.split():
                        if k.startswith("sdi"):
                            sdi = k[5:-1]
                            sdi = int(sdi)
        if sdi is None:
            abort(404)
    except:
        # raise
        abort(404)

    try:
        urlparams = urlencode({'cmd': 'text', 'type': 3, 'id': sdi, 'ofm': 3})
        url = "http://www.ncbi.nlm.nih.gov/Structure/ibis/ibis.cgi?%s" % urlparams
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
            abort(404)
    except:
        # raise
        abort(404)
    try:
        urlparams = urlencode({'type': 3, 'id': sdi})
        url = "http://www.ncbi.nlm.nih.gov/Structure/ibis/ibis.cgi?%s" % urlparams
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
        # raise
        abort(404)

    params['interac'] = interac
    params['other_sdis'] = other_sdis
    params['sdi'] = sdi
    return jsonify({'response': params})


if __name__ == "__main__":
    app = create_app()
    app.run(debug=True)

