from app.forms import *
from app.modules import *
from flask import render_template, Blueprint, request, current_app
from collections import defaultdict

interactome = Blueprint('pages', __name__)

# Routes #######

# @interactome.route("/new")
# def test_page():
#     return "TEST"


@interactome.route("/")
@interactome.route("/go_protein/<p0>/")
@interactome.route("/go_pair/<p1>/<p2>/")
@interactome.route("/go_pdb/<pdb>/")
def home_page(p0=None, p1=None, p2=None, pdb=None):

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
    return render_template('pages/home.html', species=names, error=error, p0=p0, p1=p1, p2=p2, pdb=pdb)


# @interactome.route('/')
# def home_page():
#     return render_template('pages/home.html')



# @interactome.route('/example_')
# def example_page():
#     return render_template('pages/about.html')


@interactome.route('/about')
def about_page():
    return render_template('pages/about.html')


@interactome.route("/interactome/<int:taxa>")
def interactome_page(taxa):
    return render_template('interactome.html', taxa=taxa)
    # return render_template('cy.html', taxa=taxa)


@interactome.route("/protein/<int:taxa>/<p>")
def protein_page(taxa, p):
    return render_template('protein.html', taxa=taxa, p=p)


@interactome.route("/interaction/<int:taxa>/<p1>/<p2>")
def interaction_page(taxa, p1, p2):
    return render_template('interaction.html', taxa=taxa, p1=p1, p2=p2)


@interactome.route("/interface/<int:taxa>/<p1>/<p2>/<iface>")
def interface_page(taxa, p1, p2, iface):
    return render_template('interface.html', taxa=taxa, p1=p1, p2=p2, interface=iface)


@interactome.route("/structure/<pdb>")
def structure_page(pdb):
    return render_template('structure.html', pdb=pdb)


@interactome.route("/template/<pdb>/<a>/<b>")
def template_page(pdb, a, b):
    return render_template('template.html', pdb=pdb, a=a, b=b)


###############################
#
# API
#
@interactome.route("/api/get/interactions/<int:taxa>")
def api_get_interactions(taxa):
    return jsonify({})


@interactome.route("/api/get/interactome_network/<int:taxa>")
def api_get_interactome_network(taxa):
    return jsonify({})


@interactome.route("/api/get/structure/<pdb>")
def api_get_structure(pdb):
    return jsonify({})


@interactome.route("/api/get/interaface_template/<pdb>/<a>/<b>")
def api_get_interface_template(pdb, a, b):
    return jsonify({})


# @interactome.route('/login')
# def login():
#     form = LoginForm(request.form)
#     return render_template('forms/login.html', form=form)


# @interactome.route('/register')
# def register():
#     form = RegisterForm(request.form)
#     return render_template('forms/register.html', form=form)


# @interactome.route('/forgot')
# def forgot():
#     form = ForgotForm(request.form)
#     return render_template('forms/forgot.html', form=form)
