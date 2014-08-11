import networkx as nx
import fileinput
#from pylab import *
#from scipy import *
from collections import defaultdict

###########Create 1st graph######################
G= nx.Graph()
for line in fileinput.input(["C:\Users\ragersl\Desktop\matches_human_05_20_stringent.tab"]):
 	if not fileinput.isfirstline():
		p1= line[0]
		p2= line[1]
		G.add_node(p1)
		G.add_node(p2)
		G.add_edge(p1, p2, template=line[2]) 
f.close()

##############Power Law Fitting##################
degrees_G = degree_histogram(G)
clust_G=values(clustering(G))
#clust_law= powerlaw.Fit(degrees)

#############Create 2nd graph####################
H= nx.Graph()
for line in fileinput.input(["C:\Users\ragersl\Desktop\_________.tab"]):
 	if not fileinput.isfirstline():
		p1= line[0]
		p2= line[1]
		H.add_node(p1)
		H.add_node(p2)
		H.add_edge(p1, p2, template=line[2]) 

degs_H=degree_histogram(H)
clust_H=values(clustering(H))
 
###################Comparison######################
g_edges=edges(G)
h_edges=edges(H)
union=set()
intersection=set()
for edge in g_edges:
    union.add(edge)
for edge in h_edges:
    union.add(edge)
    if edge in g_edges:
        intersection.add(edge)
jaccard=len(intersection)/len(union)

clust_array_G=[]
print "statistical significance between clustering distributions: " + scipy.stats.ttest_ind(clust_H, clust_G)[1]

##################Writing###########################

out= open("networkAnalysis", "w")
out.write("FIRST NETWORK: /n")
out.write("Number of nodes: " + G.number_of_nodes() + "\n")
out.write("Number of nodes: " + G.number_of_edges()+ "\n")
out.write("Diameter: " + diameter(G)+ "\n")
out.write("Global Clustering Coefficient: " + average_clustering(G)+ "\n")
out.write("Assortivity Coefficient (should be <0): " + degree_pearson_correlation_coefficient(G)+ "\n")
out.write("\n\n")
out.write("Jaccard Coefficient: " + jaccard + "\n")
out.close()