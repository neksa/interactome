from flask import Flask, url_for, jsonify
from flask import request, render_template
# from subprocess import call
from collections import namedtuple
from itertools import islice, ifilter

app = Flask(__name__)


Match = namedtuple('Match',
        "queryA queryB template query_type template_type " +
        "n1 m1 n2 m2 " +
        "score_template_full score_template score " +
        "scaled_score_template_full scaled_score_template scaled_score " +
        "identicalA positiveA aln_lenA bs_lenA bs_coveredA bs_alignedA bs_identicalA bs_positiveA bs_contactsA bs_BLOSUMA bs_score1A " +
        "identicalB positiveB aln_lenB bs_lenB bs_coveredB bs_alignedB bs_identicalB bs_positiveB bs_contactsB bs_BLOSUMB bs_score1B " +
        "siteA siteB")


# @app.route("/pp/", methods=['GET', 'POST'])
# def check():
#     return request.args.get('P', '?') + "-" + request.args.get('Q', '?')

"""
template resn, tpl resi, number of contacts, query resn, query resi; ....
Y,9,5,Y,193;Q,10,14,Q,194;E,78,18,E,262;R,3,36,R,187;R,74,10,R,258;I,5,10,T,189 T,46,7,P,78;Y,30,19,Y,62;E,37,12,E,69;R,49,12,R,81;A,48,9,S,80;V,33,9,I,65;F,34,25,F,66
"""


@app.route("/ppi/")
def get_list_of_protein_interactions():
    protein_pairs = [[m.queryA, m.queryB] for m in islice(iter_matches(), 0, 10)]
    return jsonify(pairs=protein_pairs)

@app.route("/ppi/<p1>/<p2>/") #if accessed without trailing slash, will redirect to page w/ slash. 
def get_interfaces(p1, p2):
    interfaces = ["{} {} {} {} {} {} {} {} {} {}".format(m.queryA, m.queryB, m.template, m.template_type, m.score_template_full, m.score_template, m.score, m.scaled_score, m.aln_lenA, m.aln_lenB) for m in ifilter(lambda x: filter_protein(x, p1, p2), iter_matches())]
    # return str(len(interfaces))
    return jsonify(interfaces=interfaces)

@app.route("/ppi/<p1>/<p2>/<interfaceID>/")
def get_interface_information(p1, p2, interfaceID): 
    interfaces = [m.__dict__ for m in ifilter(lambda x: filter_protein(x, p1, p2), iter_matches())]
    return jsonify(interfaces= interfaces)

@app.route("/ppi/<p1>/<p2>/<int:interfaceID>/alignment")
def get_interface_alignment(p1, p2, interfaceID): 
    m= list(ifilter(lambda x: filter_protein(x, p1, p2), iter_matches()))[interfaceID]
    get_residue = lambda s: {'tresn':s[0], 'tresindex':int(s[1]), 'num_contacts': int(s[2]), 'qresn':s[3], 'qresindex':int(s[4])}
    siteA= [get_residue(s.split(',')) for s in m.siteA.split(';')]
    siteB= [get_residue(s.split(',')) for s in m.siteB.split(';')]
    return jsonify(siteA=siteA, siteB=siteB)

@app.route("/")
def main_page():
    return render_template('home.html')
    # return "<form action='" + \
    #     url_for('check') + \
    #     "' method=GET>P:<input name=P type=text><p>Q:<input name=Q type=text><input type='submit' value='Check'></input></form>"

@app.route("/interface/<p1>/<p2>/")
def interface_page(p1, p2):
    return render_template('interface.html', p1=p1, p2=p2)

def iter_matches():
    fname = 'matches_human_05_20_stringent.tab'
    with open(fname) as f:
        for line in islice(f, 1, None):
            fields = line.strip().split("\t")
            m = Match._make(fields)
            yield m


def filter_protein(match, p1=None, p2=None):
    if p1 or p2:
        if match.queryA == p1 and match.queryB == p2 or match.queryA == p2 and match.queryB == p1:
            return True
    else:
        return True
    return False



if __name__ == "__main__":
    app.debug = True
    app.run()