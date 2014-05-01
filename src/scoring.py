"""
    Scoring of templates and interfaces
"""

class Scoring:
    def score(self, template, alignments, potential, profiles=None):
        """
            template is a pdb structure with two chains always named A and B (with their real names also given)
            alignment contains sequence substitutions for the ATOM sequnces of the structure, chains A and B (with protein names also given)

            the interface is detected here as distance between heavy atoms in the template
            the interface parameters are calculated in place, e.g. interface size, number of contacts, aligned residues
            all the alignment parameters like identity, similarity, gaps are recalculated here as well

            Solvent accessibility can also be calculated here

            Conservation profiles for the protins have to be provided separately
            Distance and orientation dependent potential also have to be provided 
        """

        # template is in PDB format
        # alignments is a map
        # potential is a function (distance, orientation, A_i, A_j)

        pass

