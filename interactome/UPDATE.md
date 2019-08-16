
** How to update precalculated data **

1. Fetch new PDB structures with
    update_PDB.sh

2. Activate python environment
  source env/bin/activate

  cd pipelines
  python -m interactome.workflows.analysis
  python -m interactome.workflows.temp2

3. Check for results in:
  /panfs/pan1.be-md.ncbi.nlm.nih.gov/interactomes/pipeline/Interactome/Workflow/Structures
  and
  /panfs/pan1.be-md.ncbi.nlm.nih.gov/interactomes/pipeline/Interactome/Workflow/Interfaces

