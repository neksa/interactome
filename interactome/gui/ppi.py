from flask import Flask, url_for, jsonify
from flask import request, render_template
# from subprocess import call
from collections import namedtuple
from itertools import islice

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


@app.route("/ppi/")
def get_list_of_protein_interactions():
    protein_pairs = [[m.queryA, m.queryB] for m in islice(iter_matches(), 0, 5)]
    return jsonify(pairs=protein_pairs)


@app.route("/")
def main_page():
    return render_template('home.html')
    # return "<form action='" + \
    #     url_for('check') + \
    #     "' method=GET>P:<input name=P type=text><p>Q:<input name=Q type=text><input type='submit' value='Check'></input></form>"

def iter_matches():
    fname = 'Interactome\matches_human_05_20_stringent.tab'
    with open(fname) as f:
        for line in islice(f, 1, None):
            fields = line.strip().split("\t")
            m = Match._make(fields)
            yield m


def filter_protein(match, p1=None, p2=None):
    if p1 or p2:
        if match.queryA == p1 or match.queryA == p2 or match.queryB == p1 or match.queryB == p2:
            return True
    else:
        return True
    return False



if __name__ == "__main__":
    app.debug = True
    app.run()