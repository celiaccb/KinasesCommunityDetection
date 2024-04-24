
## 0 - Package Initialization

#Packages necessary to construct a kinase network from the user uploaded file
import pandas as pd
import networkx as nx
from networkx import *
import csv
import itertools
from itertools import chain
import numpy as np
from scipy.stats import norm
import statistics
from statistics import *
import operator
import functools
#Visualizations
import plotly.express as px
import plotly.graph_objects as go 
from pyvis.network import Network
#Community detection
import community as community_louvain
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from sknetwork.clustering import Louvain
#import random
import collections

## 1 - Functions

### 1.1 - Extracting information from a z-scores dataset

# input: dataset must be a csv files in which the first column is the kinase interactions (pairs), and each following 
# column the z-scores of each of the kinase interactions in one of the treatments
# output: treatments are the headers of each of the z-scores columns, used to name visulizations returned from the analysis
# negative_edges is a list of lists, each list corresponding to a treatment in the dataset, it's contents being the 
# kinase interactions in said treatment with a z_score <0 and their z_scores,in the shape of tuples, i.e. a network
# negative_nodes is a list of lists, each list corresponding to the kinases of one network
def dataset_info(dataset):
    # Read dataset and drop rows with no z-score information
    networkfile = pd.read_csv(dataset)
    networkfile = networkfile.dropna()
    
    # Return the headers of the z_scores columns, i.e. the treatments names
    headers = list(networkfile.columns.values)
    treatments = headers[1:]
    
    # Convert first column (kinase interactions/pairs) of the document to a list
    edges_list= networkfile[networkfile.columns[0]].to_list()
    # Convert the list into a list of tuples, each tuple being one kinase pair/edge (since this is the format needed to 
    # add edges onto pyvis.network (network visualization package))
    edges_list2 = [tuple(x.split('.',1)) for x in edges_list]
    
    #For each treatment(column) on the file, store kinase interactions/pairs with a z-score <0, and their z-scores
    negative_edges = [[(*edges_list2[n],+ networkfile[networkfile.columns[na]].to_list()[n]) for n in range(len(edges_list2)) if networkfile[networkfile.columns[na]].to_list()[n] < 0 ] for na in range(1,len(headers))]
    
    #Kinases in each of the negative_edges lists
    negative_nodes = []
    for selected_treatment in negative_edges:
        #Keep only the kinases names (i.e. remove the z-scores) from the list of edges
        nodes_bf = [ (a,b) for a,b,c in selected_treatment]
        #Flatten list of edges tuples to list of nodes
        nodes = list(itertools.chain(*nodes_bf))
        #Remove duplicates from the list
        nodes_list = list(dict.fromkeys(nodes))
        negative_nodes.append(nodes_list)
        
    return treatments,negative_edges, negative_nodes


### 1.2 - Community detection

# input: edges_list is a list corresponding to a treatment in the dataset, it's contents being the kinase interactions
# in said treatment with a z_score <0 and their z_scores,in the shape of tuples, i.e. a network 
# output: A_asso is the converged association matrix, n_iter the number of iterations required for convergence and 
# labels are the converged community labels 
def converge_asso(edges_list):
    
    #Initialize a Networkx instance
    negativenetwork = nx.Graph()
    #Add the kinase interactions & z-scores of treatment w to the network
    negativenetwork.add_weighted_edges_from(edges_list)
    #A_0 is the adjacency matrix of the matrix (network instance) stored as a numpy array
    A_0 = to_numpy_array(negativenetwork)
    
    n_runs = 100
    A_runs = np.zeros((n_runs, len(A_0)))
    #Run the community detection algorithm on the network A_0 an n_runs amount of times
    for i in range(n_runs):
        #Louvain algorithm with Potts null network
        louvain = Louvain(1, 'potts', shuffle_nodes = 1)
        #Communities identified by number (Community 1, Community 2...)
        #Array is the length of the number of kinases in the network
        #For each kinase, it's community assignement is indicated by number
        labels = louvain.fit_transform(A_0)
        #Store the community partitions of each run
        A_runs[i,:] = labels

    #X = Number of kinases in the network
    #Create & store 'X' arrays of 'X' zeros each
    A_asso = np.zeros((len(A_0), len(A_0)))
    #After loop, each array corresponds to one kinase, and each position within an array indicates in how many runs that kinase is
    #assigned to the same community as each of the kinases in the network 
    #Loop over each run
    for i in range(n_runs):
        #Loop over each kinase/node in network (identified by number)
        for j in range(len(A_0)):
            #Loop over each kinase/node in network (identified by number)
            for l in range(len(A_0)):
            #If j and l are assigned to the same community in run i, add '1' to position l of array j   
                A_asso[j,l] = A_asso[j,l] + (A_runs[i,j] == A_runs[i,l])
    
    #Divide each element in each of the arrays in the converged association matrix(A_asso) by n_runs to get
    #the percentage (in 0.1(=10%) format) of runs in which each kinase-kinase/node-node combination is assigned to the same community
    A_asso = np.true_divide(A_asso, n_runs)
    
    n_iter = 0
    #until all the entries in A_asso are 0(not grouped together) or 1(grouped together) (all runs are the same)
    while(len(np.unique(A_asso)) != 2):
        n_runs = 100
        #Create an array for each run, the length of each array being the number of kinases in the network
        A_runs = np.zeros((n_runs, len(A_0)))
        #Repeat community detection on co-classification matrix (the weight of the edges is the percentage of runs in which two kinases are grouped together)
        for i in range(n_runs):
            louvain = Louvain(1, 'potts', shuffle_nodes = 1)
            labels = louvain.fit_transform(A_asso)
            A_runs[i,:] = labels
        A_asso = np.zeros((len(A_0), len(A_0)))
        for i in range(n_runs):
            for j in range(len(A_0)):
                for l in range(len(A_0)):
                    A_asso[j,l] = A_asso[j,l] + (A_runs[i,j] == A_runs[i,l])
        #Return percentage of runs in which each kinase-kinase combination is a ssigned to the same community
        A_asso = np.true_divide(A_asso, n_runs)
        n_iter = n_iter + 1 
    
    #Final community partitions (all runs are the same)
    labels = A_runs[0,:]                                           
    return A_asso,n_iter,labels

### 1.3 - Storing the contents of the community containing a selected kinase

# input: edges_list is a list corresponding to a treatment in the dataset, it's contents being the kinase interactions
#in said treatment with a z_score <0 and their z_scores,in the shape of tuples, i.e. a network
# target is the main targeted kinase of a treatment/inhibitor (i.e. the selected kinase), labels is 
# output: community is the community in the network that contains the targeted kinase
def community_contents(edges_list, target, labels):
    
    #Initialize a Networkx instance
    negativenetwork = nx.Graph()
    #Add the kinase interactions & z-scores of treatment w to the network
    negativenetwork.add_weighted_edges_from(edges_list)
    #Store the individual kinases in the network as a numpy array
    node_0 = list(negativenetwork.nodes())
    node_0_arr = np.asarray(node_0)
    #Convert network to a numpy array
    A_0 = to_numpy_array(negativenetwork)
        
    #Find the indice of the kinase of interest within the kinases array
    target_ind = np.argwhere(node_0_arr == target)
    #Find the community assignement of the target (number), and return the indices of all the kinases with that number (i.e. same community)
    com_cont_labels = np.argwhere(labels == int(labels[target_ind]))
    #Use said indices to find the name of the community kinases in the nodes array and store them
    community = node_0_arr[com_cont_labels]
    #Convert to list of kinases
    community = list(itertools.chain(*community))
    return community


### 1.4 - Reproducibility check

#Run converge_asso an n_checks number of times and check if the contents of the community of interest differ
def rep_check(n_checks, edges_list, target):
    community = ['a'] * n_checks
    for n in range(n_checks):
        #Store community partitions of this check
        A_checks = converge_asso(edges_list)[2]
        #Store contents of community of interest
        community[n] = community_contents(edges_list,target,labels=A_checks)
    #Check if there are any kinases that are not present in all the community versions
    all_com = list(itertools.chain(*community))
    kin_rem = [x for x in set(all_com) if all_com.count(x)==n_checks]
    return kin_rem

### 1.5 - Network visualization

# input: edges_list is a list corresponding to a treatment in the dataset, it's contents being the kinase interactions
# in said treatment with a z_score <0 and their z_scores,in the shape of tuples, i.e. a network
# kinase_groups is a list of the groups in which we want to classify the kinases in the network
# kg_colors is a list of the colors in which each kinase groups should appear in the visualization
# treatment is the header of the column from which the z-scores for edges_list was extracted (i.e. name of the treatment)
# output: A visualization of the network constructed from the edges list, in HTML format
def network_visualization(edges_list, kinases, kg_colors, treatment):    
    #Initialize a pyvis.network instance
    networkgraph = Network(height='1000px',width='1000px', heading='')
    
    #Add each community of nodes with their assigned color
    for n in range(len(kinases)):
        for x in kinases[n]:
            networkgraph.add_node(x, value=15000000000, title=x, label=x, color=kg_colors[n], mass=23)
    
    #Add the edges
    networkgraph.add_edges(edges_list)
    
    #Make empty dictionary in which to store the degree of each kinase
    degree_graphs = {}
    #Save list of nodes to extract information such as kinases degree centrality
    kinases_map = networkgraph.get_adj_list()

    for node in networkgraph.nodes:
        #Add degree to each kinase in the empty dictionary
        degree_graphs[node["title"]] = len(kinases_map[node["id"]]) 
        #Add box showing degree + to which nodes each node is connected, which appears when the user hovers over a node on the network graph
        node["title"] += " connected to "+ str(len(kinases_map[node["id"]])) +" kinases:<br>" + "<br>".join(kinases_map[node["id"]]) 
        kinase = len(kinases_map[node["id"]])
        #Change the nodes to be bigger the higher their degree is
        if  0 <=  kinase <= 19 :
            node["value"] = 100000

        if  20 <=  kinase <= 29 :
            node["value"] = 1000000000

        if  30 <=  kinase <= 39 :
            node["value"] = 2000000000

        if  40 <=  kinase <= 49:
            node["value"] = 10000000000
    
        if  50 <=  kinase <= 59:
            node["value"] = 20000000000

        if 60 <=  kinase <= 69:
            node["value"] = 30000000000
        
        if  70 <=  kinase:
            node["value"] = 40000000000
    
    #Settings regarding the appearance of the network
    networkgraph.set_options("""
    var options = {
    "nodes": {
    "borderWidth": 1.5,
        "font": {
        "size": 80,
        "face": "verdana"
        }
      },
      "edges": {
        "color": {
          "color": "rgba(192,187,223,1)",
          "inherit": true
        },
        "font": {
          "strokeWidth": 6
        },
        "hoverWidth": 3.2,
        "smooth": false
      },
      "physics": {
        "barnesHut": {
          "springLength": 1040
        },
        "minVelocity": 0.75
      }
    }
    """)
    
    #Return visualization as HTML file
    final_network = treatment + '[-]'
    networkgraph.show('output_visualizations/network_visualization{}.html'.format(final_network))

### 1.6 - Community strength calculation & display in heatmap

# input: edges_list is a list corresponding to a treatment in the dataset, it's contents being the kinase interactions
# in said treatment with a z_score <0 and their z_scores,in the shape of tuples, i.e. a network
# community is a list of the kinases in the community of interest
# output: dictionary of kinases and their community strength
def cs_calculation(treat_int):
    #Calculate the community strength of each of the kinases in the two selected communities you are comparing
    kin_cs = ['a'] * len(treat_int)
    for n in range(len(treat_int)):
        na = treat_int[n]
        # Store kinase pairs in which both kinases are in the community of interest
        community_edges = [x for x in negative_edges[na] if x[0] in communities[na] and x[1] in communities[na]]
        # Calculate the community strength of each of the kinases in the community and store as dictionary
        kin_cs[n] = {kin:abs(round(sum([community_edges[n][2] for n in range(len(community_edges)) if community_edges[n][1]==kin or community_edges[n][0]==kin]), 2)) for kin in communities[na]}    
    
    #Combine two communities lists
    kinases = list(dict.fromkeys(communities[treat_int[0]]+communities[treat_int[1]]))
    #Get values for each of the kinases in the combined list in each of the communities
    #If kinase is not in the community then the value is None
    values = [[kin_cs[n][kin] if kin in communities[treat_int[n]] else None for kin in kinases]for n in [0,1]]
    
    #Define heatmap axis and text
    x = kinases
    y = [ treatments[na]+' community CS' for na in treat_int]
    z = [values[0],values[1]]
    
    #kinases, treatment,  = list(zip(x, y, z))
    
    return [x,y,z]
    
 ## 2 - Analysis and Results

 #Store dataset treatments names, kinase interactions (network edges) with a z-score < 0, and their z-scores (edges weight) 
[treatments,negative_edges, negative_nodes] = dataset_info('input/PROJECT_DATASET_2.csv')
#Main kinase targets of each of the treatments/inhibitors in the dataset
targets = ['MAP2K1','PIK3CA','AKT1_2','MAPK1_3']

### 2.1 - Selected communities contents for networks trametinib(MAP2K1)[-], GDC0994(MAPK1_3), GDC0941(PI3K)[-], and AZD5363(AKT)[-]

#Please note that the following block of code takes a few minutes to run

#Store contents of communities of interest in each network
communities = ['a'] * 4
for n in [1,2,3]:
    [A_asso,n_iter,labels] = converge_asso(negative_edges[n])
    communities[n] = community_contents(negative_edges[n], targets[n], labels)      
#Because we noticed that when extracting the community of interest from the trametinib(MEK)[-] network with the converge_asso
#function the contents of said community differ from run from run, for this network we run the community detection function (converge_asso)
#50 times, and keep only the nodes that are present in the community in all 50 runs of the converge_asso function
communities[0] = rep_check(50, negative_edges[0], targets[0])

#Community containing MAP2K1 in the trametinib(MEK)[-] network (column 1 , negative z-scores)
#Community containing MAPK1_3 in the GDC0994(ERK)[-] network (column 4, negative z-scores)
MEK_ERKcom = [communities[0],communities[3]]
#Separate kinases in networks trametinib(MEK)[-] and GDC0994(ERK)[-] into:
#Kinases not present in the community of interest of that network
MEK_ERKnonc = [[x for x in negative_nodes[0] if x not in communities[0]],[x for x in negative_nodes[3] if x not in communities[3]]]
#Kinases present in both communities
MEK_ERKallc = [x for x in communities[0] if x in communities[3]]
MEK_ERKallc = [MEK_ERKallc,MEK_ERKallc]
#Kinases present in one of the communities, but not the other
MEK_ERKinc = [[x for x in communities[0] if x not in communities[3]],[x for x in communities[3] if x not in communities[0]]]

#Community containing PIK3CA in the GDC0941(PI3K)[-] network (column 2, negative z-scores)
#Community containing AKT1_2 in the AZD5363(AKT)[-] network (column 3, negative z-scores)
PI3K_AKTcom = [communities[1], communities[2]]
#Separate kinases in networks GDC0941(PI3K)[-] and AZD5363(AKT)[-] into:
#Kinases not present in the community of interest of that network
PI3K_AKTnonc = [[x for x in negative_nodes[1] if x not in communities[1]],[x for x in negative_nodes[2] if x not in communities[2]]]
#Kinases present in both communities
PI3K_AKTallc = [x for x in communities[1] if x in communities[2]]
PI3K_AKTallc = [PI3K_AKTallc, PI3K_AKTallc]
#Kinases present in one of the communities, but not the other
PI3K_AKTinc = [[x for x in communities[1] if x not in communities[2]],[x for x in communities[2] if x not in communities[1]]]

nonc = [MEK_ERKnonc[0],PI3K_AKTnonc[0],PI3K_AKTnonc[1],MEK_ERKnonc[1]]
allc = [MEK_ERKallc[0],PI3K_AKTallc[0],PI3K_AKTallc[1],MEK_ERKallc[1]]
inc = [MEK_ERKinc[0],PI3K_AKTinc[0],PI3K_AKTinc[1],MEK_ERKinc[1]]
        
kinase_groups = list(zip(nonc,allc,inc))
kinase_int = list(zip(allc,inc))

#Show the community contents of the selected communities in networks trametinib(MAP2K1)[-], GDC0941(PI3K)[-],
#AZD5363(AKT)[-] and GDC0994(MAPK1_3), in that particular order
for n in range(len(communities)):
    print(treatments[n])
    print(communities[n])

 ### 2.2 - Network visualizations of networks trametinib(MAP2K1)[-], GDC0994(MAPK1_3), GDC0941(PI3K)[-], and AZD5363(AKT)[-]

for n in range(len(communities)):
    networkgraph = network_visualization(edges_list=negative_edges[n],kinases=kinase_groups[n],kg_colors=["#fcbb8b", "#857be3", "#baeeff"], treatment=treatments[n])
    
### 2.3 - Network visualization of intersection between the selected communities of trametinib(MAP2K1)[-], GDC0941(PI3K)[-] and AZD5363(AKT)[-]

#Kinases present in the selected community of network trametinib(MAP2K1)[-], but not in the selected communities of 
#networks GDC0941(PIK3CA)[-] and AZD5363(AKT1_2)[-]
#MAP2K1only = [x for x in communities[0] if x not in communities[1] and x not in communities[2]]
#Kinases present accross the selected communities of networks trametinib(MAP2K1)[-], GDC0941(PIK3CA)[-] and AZD5363(AKT1_2)[-]
intersection = [x for x in communities[0] if x in communities[1] and x in communities[2]]
#Kinases present accross the selected communities of networks trametinib(MAP2K1)[-] and GDC0941(PIK3CA)[-], but not in the
#selected community of network AZD5363(AKT)[-]
MAP2K1_PIK3CA = [x for x in communities[0] if x in communities[1] and x not in communities[2]]
#Kinases present accross the selected communities of networks AZD5363(AKT)[-] and GDC0941(PIK3CA)[-] but not trametinib(MAP2K1)[-]
#and kinases present only in the selected communities of GDC0941(PIK3CA)[-] and AZD5363(AKT1_2)[-] networks
#therest = [[x for x in communities[3] if x in communities[2] and x not in communities[0]],[x for x in communities[3] if x not in communities[2] and x not in communities[0]], [x for x in communities[2] if x not in communities[3] and x not in communities[0]]]

kinase_groups = MAP2K1_PIK3CA, intersection
kinases = list(itertools.chain(*kinase_groups))
#Store kinase interactions across the three networks in which the two kinases are in any of the selected communities
#edges = negative_edges[0] + negative_edges[1] + negative_edges[2]
#edges = [(a,b,0) for a,b,c in edges if a in kinases and b in kinases]
edges = [ [(a,b) for a,b,c in negative_edges[n]] for n in [0,1,2] ]
edges = set(edges[0]).intersection(edges[1], edges[2])
edges = [(a,b) for a,b in edges if a in kinases and b in kinases]

#Intersection of communities visualized as a network
networkgraph= network_visualization(edges_list=edges,kinases=kinase_groups,kg_colors=["#fcbb8b", "#857be3", "#baeeff","green"], treatment='intersection_MAP2K1_PIK3CA_AKT')

### 2.4 - Community strength of the kinases in each of the selected communities of networks trametinib(MAP2K1)[-], GDC0941(PI3K)[-],AZD5363(AKT)[-] and GDC0994(MAPK1_3), visualized in heatmaps

#trametinib(MAP2K1)[-] and GDC0994(MAPK1_3)
[x1,y1,z1] = cs_calculation(treat_int=[0,3])

# HEATMAP
fig = go.Figure(data=go.Heatmap(
                z=z1,
                x=x1,
                y=y1,
                hoverongaps = False))
fig.show()

#GDC0941(PI3K)[-] and AZD5363(AKT)[-]
[x2,y2,z2] = cs_calculation(treat_int=[1,2])

# HEATMAP
fig = go.Figure(data=go.Heatmap(
                z=z2,
                x=x2,
                y=y2,
                hoverongaps = False))
fig.show()

### 2.5 - Saving the edges and nodes of each community as dataframes

# Individual communities

#Save nodes and edges as dataframe to make nice visualizations in R

#nodes community strength
nodes_cs = list(zip([x1]+[x1],z1)) + list(zip([x2]+[x2],z2)) 
nodes_cs = [nodes_cs[0]] + [nodes_cs[2]] + [nodes_cs[3]] + [nodes_cs[1]]

#Filter to keep edges of interest only
negative_edges_com = [ [x for x in negative_edges[n] if (x[0] in list(itertools.chain(*kinase_int[n])) and x[1] in list(itertools.chain(*kinase_int[n])) ) ] for n in range(4)]
#Add weight to draw edges
negative_edges_pos = [ [ (*x , abs(x[2])*10 ) for x in network] for network in negative_edges_com ]

#Add community + color to nodes
nodes_inf = []
for n in range(len(kinase_int)):
    nodes_inf0 = list(zip( kinase_int[n][0], ['kinases present exclusively in this SC']*len(kinase_int[n][0]) ,  ['#96efff']* len(kinase_int[n][0]) )) 
    nodes_inf1 = list(zip(kinase_int[n][1], ['kinases present in both SCs']*len(kinase_int[n][1]) ,  ['#FFDB33']* len(kinase_int[n][1]) )) 
    nodes_inf_sum = nodes_inf0 + nodes_inf1
    nodes_df = pd.DataFrame( nodes_inf_sum, columns=['id', 'group', 'color'])
    #add community strength
    nodes_cs_df = pd.DataFrame( list(zip(nodes_cs[n][0], nodes_cs[n][1])) , columns = ['id','community_strength'] )
    nodes_df_add = pd.merge(nodes_df, nodes_cs_df, on='id')
    #save nodes as csv file
    nodes_df_add.to_csv('output_dataframes/{}_nodes.csv'.format(treatments[n]))
    #save edges as csv file
    edges_df = pd.DataFrame(negative_edges_pos[n], columns=['from','to','z-score','weight'])
    edges_df.to_csv('output_dataframes/{}_edges.csv'.format(treatments[n]))


# Intersection of trametinib, GDC0941 and AZD5363

#Calculate average edge weights across the three networks
edges_MAP2K1_PIK3CA_AKT = [ { (edge[0],edge[1]):abs(edge[2]) for edge in negative_edges[n]} for n in [0,1,2]  ]   
edges_MAP2K1_PIK3CA_AKT_sum = dict( functools.reduce(operator.add, map(collections.Counter, xx))  )
edges_MAP2K1_PIK3CA_AKT_ave = [ (*key, (result[key]/3)) for key in edges ]

#Save nodes and edges as dataframe to make nice network visualizations in R
edges_df = pd.DataFrame(edges_MAP2K1_PIK3CA_AKT_ave, columns=['from','to','weight'])
edges_df.to_csv('output_dataframes/MAP2K1_PIK3CA_AKT_edges.csv')

intersection_x = list(zip(intersection, ['MAP2K1_PIK3CA_AKT'] * len(intersection), ['#96efff']* len(intersection) ) )
MAP2K1_PIK3CA_x = list(zip(MAP2K1_PIK3CA, ['MAP2K1_PIK3CA_notAKT']*len(MAP2K1_PIK3CA), ['#FFDB33']* len(MAP2K1_PIK3CA) ) )
nodes = intersection_x + MAP2K1_PIK3CA_x
nodes_df = pd.DataFrame( nodes, columns=['id', 'group', 'color'])

nodes_df.to_csv('output_dataframes/MAP2K1_PIK3CA_AKT_nodes.csv')





