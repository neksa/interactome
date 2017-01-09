from flask import Flask, url_for, jsonify
from flask import request, render_template, redirect
from sqlalchemy import create_engine, MetaData, Table
# from subprocess import call
from collections import namedtuple
from itertools import islice, ifilter

app = Flask(__name__)

engine = create_engine('mysql://localhost/agoncear', convert_unicode=True, echo=True)
metadata = MetaData(bind=engine)

species = Table('species', metadata, autoload=True)
pairs = Table('pairs', metadata, autoload=True)
proteins = Table('proteins', metadata, autoload=True)
idmapping = Table('idmapping', metadata, autoload=True)

# @app.route("/pp/", methods=['GET', 'POST'])
# def check():
#     return request.args.get('P', '?') + "-" + request.args.get('Q', '?')

"""
template resn, tpl resi, number of contacts, query resn, query resi; ....
Y,9,5,Y,193;Q,10,14,Q,194;E,78,18,E,262;R,3,36,R,187;R,74,10,R,258;I,5,10,T,189 T,46,7,P,78;Y,30,19,Y,62;E,37,12,E,69;R,49,12,R,81;A,48,9,S,80;V,33,9,I,65;F,34,25,F,66
"""


@app.route("/test/")
def test():
    return render_template('cy.html')


@app.route("/")
@app.route("/go_protein/<p0>/")
@app.route("/go_pair/<p1>/<p2>/")
@app.route("/go_pdb/<pdb>/")
def intro_page(p0=None, p1=None, p2=None, pdb=None):

    def is_protein(p):
        return len(p) == 5

    def is_complex(pdb):
        return len(pdb) == 4

    error = None
    if p0 is not None:
            if is_protein(p0):
                return redirect("/protein/{}/".format(p0))
            else:
                error = "Not a valid protein identifier: {}".format(p0)

    if p1 is not None and p2 is not None:
        if is_protein(p1) and is_protein(p2):
            return redirect("/interaction/{}/{}".format(p1, p2))
        else:
            error = "One of the two proteins is not a valid protein identifier: {}, {}".format(p1, p2)

    if pdb is not None:
        if is_complex(pdb):
            return redirect("/structure/{}/".format(pdb))
        else:
            error = "PDB ID not valid or not a protein complex: {}".format(pdb)

    # get list of species
    rows = species.select().execute()
    names = {}
    for r in rows:
        names[r.name] = r.tax_id
    return render_template('intro.html', species=names, error=error, p0=p0, p1=p1, p2=p2, pdb=pdb)


@app.route("/interactome/<int:taxa>/")
def page_interactome(taxa):
    return render_template('interactome.html', taxa=taxa)
    # return render_template('cy.html', taxa=taxa)


@app.route("/protein/<int:taxa>/<p>/")
def page_protein(taxa, p):
    return render_template('protein.html', taxa=taxa, p=p)


@app.route("/interaction/<int:taxa>/<p1>/<p2>/")
def page_interaction(taxa, p1, p2):
    return render_template('interaction.html', taxa=taxa, p1=p1, p2=p2)


@app.route("/interface/<int:taxa>/<p1>/<p2>/<iface>/")
def page_interface(taxa, p1, p2, iface):
    return render_template('interface.html', taxa=taxa, p1=p1, p2=p2, interface=iface)


@app.route("/structure/<pdb>/")
def page_structure(pdb):
    return render_template('structure.html', pdb=pdb)


@app.route("/template/<tpl>/")
def page_template(tpl):
    return render_template('template.html', tpl=tpl)


# @app.route("/cytoscape/")
# def page_cytoscape():
#     return render_template('cytoscape.html')


####################################

@app.route("/ppi/<int:taxa>/")
def ajax_list_interactions(taxa):
    rows = pairs.select().where(
        pairs.c.tax_id == taxa).where(
        pairs.c.bs_positive_a > 20).where(
        pairs.c.bs_positive_b > 20).limit(1000).execute().fetchall()
    protein_pairs = [[m.uniprot_ac_a, m.uniprot_ac_b] for m in rows]
    return jsonify(pairs=protein_pairs)


@app.route("/ppi/cy/<int:taxa>/")
def ajax_list_interactions_cy(taxa):
    rows = pairs.select().where(
        pairs.c.tax_id == taxa).where(
        pairs.c.bs_positive_a > 20).where(
        pairs.c.bs_positive_b > 20).limit(10000).execute().fetchall()
    protein_pairs = [(m.pair_id, m.uniprot_ac_a, m.uniprot_ac_b) for m in rows]
    cy = []
    cy.extend([{'group': 'nodes', 'data': {'id': n, 'name': n}} for i, n, m in protein_pairs])
    cy.extend([{'group': 'nodes', 'data': {'id': m, 'name': m}} for i, n, m in protein_pairs])
    cy.extend([{'group': 'edges', 'data': {'id': str(i), 'source': n, 'target': m}} for i, n, m in protein_pairs])
    # cy = {'nodes': nodes, 'edges': edges}
    # return jsonify(graph=cy)
    return jsonify(graph=cy)


@app.route("/ppi/<int:taxa>/<p1>/<p2>/")  # if accessed without trailing slash, will redirect to page w/ slash.
def ajax_get_interfaces(taxa, p1, p2):
    rows = pairs.select().where(
        pairs.c.tax_id == taxa).where(
        pairs.c.uniprot_ac_a == p1).where(
        pairs.c.uniprot_ac_b == p2).execute().fetchall()
    interfaces = [dict(zip(i.keys(), i.values())) for i in rows]
    return jsonify(interfaces=interfaces)


@app.route("/ppi/<int:taxa>/<p1>/<p2>/<int:pair_id>/")
def ajax_get_interface_information(taxa, p1, p2, pair_id):
    row = pairs.select().where(
        pairs.c.tax_id == taxa).where(
        pairs.c.uniprot_ac_a == p1).where(
        pairs.c.uniprot_ac_b == p2).where(
        pairs.c.pair_id == pair_id).execute().fetchone()
    properties = dict(zip(row.keys(), row.values()))
    return jsonify(properties=properties)


@app.route("/ppi/<int:taxa>/<p1>/<p2>/<int:pair_id>/alignment")
def ajax_get_interface_alignment(taxa, p1, p2, pair_id):
    row = pairs.select().where(
        pairs.c.tax_id == taxa).where(
        pairs.c.uniprot_ac_a == p1).where(
        pairs.c.uniprot_ac_b == p2).where(
        pairs.c.pair_id == pair_id).execute().fetchone()
    properties = dict(zip(row.keys(), row.values()))
    s_a = properties['site_a']
    s_b = properties['site_b']
    print s_a
    print s_b
    a1 = ""
    a2 = ""
    for r in s_a.split(';'):
        i = r.split(',')
        a1 += i[0]
        a2 += i[3]
    a = a1 + "<br>" + a2

    b1 = ""
    b2 = ""
    for r in s_b.split(';'):
        i = r.split(',')
        b1 += i[0]
        b2 += i[3]
    b = b1 + "<br>" + b2
    print a
    print b

    return jsonify(site_a=a, site_b=b)
    # m = list(ifilter(lambda x: filter_protein(x, p1, p2), iter_matches()))[interfaceID]
    # get_residue = lambda s: {'tresn': s[0], 'tresindex': int(s[1]), 'num_contacts': int(s[2]), 'qresn': s[3], 'qresindex': int(s[4])}
    # siteA = [get_residue(s.split(',')) for s in m.siteA.split(';')]
    # siteB = [get_residue(s.split(',')) for s in m.siteB.split(';')]
    # return jsonify(siteA=siteA, siteB=siteB)


# def iter_matches():
#     fname = 'matches_human_05_20_stringent.tab'
#     with open(fname) as f:
#         for line in islice(f, 1, None):
#             fields = line.strip().split("\t")
#             m = Match._make(fields)
#             yield m


# def filter_protein(match, p1=None, p2=None):
#     if p1 or p2:
#         if match.queryA == p1 and match.queryB == p2 or match.queryA == p2 and match.queryB == p1:
#             return True
#     else:
#         return True
#     return False


if __name__ == "__main__":
    app.debug = True
    app.run()
