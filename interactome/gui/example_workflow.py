import json
from urllib.request import urlopen

prefix = "https://www.ncbi.nlm.nih.gov/projects/Interactomes/interactomes/rest/v1"
API_similar_compound = "{}/compounds/similar/compound/{}"
API_structures_compound = "{}/structures/compound/{}"
API_obs_prot_compound_bs = "{}/sites/obs/compound/pdb/{}/{}"

combined_results = []

# Find similar compounds
query_cid = 11433190  # Nultin-3a
URI1 = API_similar_compound.format(prefix, query_cid)
with urlopen(URI1) as f1:
    result_cids = json.loads(f1.read().decode('utf-8'))
    for cid in result_cids['response']:

        # Find structures for each compound
        URI2 = API_structures_compound.format(prefix, cid)
        with urlopen(URI2) as f2:
            result_struct = json.loads(f2.read().decode('utf-8'))
            for pdb, chain in result_struct['response']:

                # Find small molecules and binding sites for each structure
                URI3 = API_obs_prot_compound_bs.format(prefix, pdb, chain)
                with urlopen(URI3) as f3:
                    results_sites = json.loads(f3.read().decode('utf-8'))
                    combined_results.extend(results_sites['response'])

with open("results.json", 'w') as output:
    output.write(json.dumps(combined_results))
