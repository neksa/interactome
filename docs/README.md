PDB templates:

195829 interfaces
84275 if we only count the ones with at least one nongenerated chain
20690 out of them contains one generated chain


"dimers" pdb+chain+chain in PDB with well-defined interfaces (5res - 5res with <4A contact): 60268
pdb codes with well defined interfaces: 23913




iterator on atom entries (they have correct numbering)
 chain+model
 resn+resi
 atomn
 xyz




Data formats:

   I. */atomic_contacts_5A.tab
  II. */interface_residue_atoms_5A.tab
  
 III. */interface_properties_4A.tab
  IV. */interface_residues_contacts_4A.tab
   V. distance_stats_4A.tab distance_stats_reference_4A.tab
  VI. potential_4A.index

 VII. ?? other scoring data, like BLOSUM62
VIII. scoring of interface_matches_4A.tab
  IX. network, each interface (MITAB 2.7)
   X. network, gene-gene, one interaction (MITAB 2.7)



Info: Number of templates (PDB+chain+chain) =  155311

A problematic case: 1MYP chain C, completely wrong matching between seqres and atomres. Why?
Another case: 2PYO chains C and D
3O4O B
3QUM B
4JY4 A, B
1AOI D
1D4V A
1R5V B
3D3W B, B_1
3BX7 A


the interface residue in template hit (SEQRES) 1AOI D 10 S does not match what's in the structure (ATOMRES): N 10 10
alignment, pos = 4
TRKESYAIYVYKVLKQVHPDTG---ISSKAMSIMNSFVNDVFERIAGEASRLAHYNKRSTITSREIQTAVR

the interface residue in template hit (SEQRES) 1AOI D 15 V does not match what's in the structure (ATOMRES): T 15 15
alignment, pos = 9
TRKESYAIYVYKVLKQVHPDTG---ISSKAMSIMNSFVNDVFERIAGEASRLAHYNKRSTITSREIQTAVR

the interface residue in template hit (SEQRES) 1AOI D 13 I does not match what's in the structure (ATOMRES): G 13 13
alignment, pos = 7
TRKESYAIYVYKVLKQVHPDTG---ISSKAMSIMNSFVNDVFERIAGEASRLAHYNKRSTITSREIQTAVR


