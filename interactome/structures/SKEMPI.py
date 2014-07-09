
class SKEMPI:
    def getDimers(self, fname="../SKEMPI/dimers.csv"):
        dimers = defaultdict(list)
        with open(fname) as f:
            for line in f:
                pdb, chain1, chain2 = line.strip().split("_", 2)
                dimers[pdb] = (chain1, chain2)
        return dimers
