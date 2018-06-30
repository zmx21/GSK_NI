import louvain
import igraph as ig
import sys
#fileName input argument, should be path to edge file (3 cols: Node1, Node2, Weight)
def run_louvain(fileName):
	#Construct igraph
	g = ig.Graph.Read_Ncol(fileName,names=True,weights=True,directed=False)
	#Find clusters, using louvain. Pass in weights that's same order as edges. 
	partition = louvain.find_partition(g,louvain.ModularityVertexPartition,weights=g.es["weight"])
	#print(g.vs.indices)
	#print(vars(partition))
	#Store vertices info
	vertices = g.vs
	#Get clusters
	membershipList = partition.membership
	#Get Names for vertices that matches order of membership
	verticeNames = []
	for i in range(0,len(vertices)):
		verticeNames.append(vertices[i]["name"])
	#print(membershipList)
	#print(verticeNames)
	return membershipList,verticeNames	

#Write to output file, first column is the node name. Second column is the cluster number. 
#fo = open(fileName.split('.')[0]+".out", "w");
#for i in range(0,len(vertices)):
	#print(str(vertices[i]))
	#fo.write( vertices[i]["name"]+" "+str(membershipList[i])+"\n" );
#fo.close();
