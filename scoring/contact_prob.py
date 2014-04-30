#
# Contact probabilities
#


with open("aa-distances.tab", 'r') as f:
    for line in f:
        aa, d12, dCA12 = line.strip().split(3)
        
