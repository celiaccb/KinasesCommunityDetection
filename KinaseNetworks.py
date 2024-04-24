
0 - Package Initialization------------------------------------------------------

##Python version used: 3.9.12
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
#from pyvis.network import Network
import matplotlib as mpl
import matplotlib.colors as mcolors
#Community detection
import community as community_louvain
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from sknetwork.clustering import Louvain
#Data manipulation
import collections
from collections import ChainMap
from collections import OrderedDict
import BioVenn
from BioVenn import draw_venn
import re
from itertools import chain,combinations
from funcy import join_with

#---------------------------------------------------------------------------------

#1 - Functions--------------------------------------------------------------------

#1.1 - Extracting information from a z-scores dataset

# input: dataset must be a csv files in which the first column is the kinase interactions (pairs), and each following 
# column the z-scores of each of the kinase interactions in one of the treatments
# output: treatments are the headers of each of the z-scores columns, used to name visulizations returned from the analysis
# negative_edges is a list of lists, each list corresponding to a treatment in the dataset, it's contents being the 
# kinase interactions in said treatment with a z_score <0 and their z_scores,in the shape of tuples, i.e. a network
# negative_nodes is a list of lists, each list corresponding to the kinases of one network
def dataset_info(dataset, rem_vs_control='no', sign = 'neg'):
    # Read dataset and drop rows with no z-score information
    networkfile = pd.read_csv(dataset)
    networkfile = networkfile.dropna()
    
    # Return the headers of the z_scores columns, i.e. the treatments names
    headers = list(networkfile.columns.values)[1:]
    treatments = [x.split(".")[-1] for x in headers]
    if rem_vs_control == 'yes':
        treatments = [x.split("_")[0] for x in treatments]
    
    # Convert first column (kinase interactions/pairs) of the document to a list
    edges_list= networkfile[networkfile.columns[0]].to_list()
    # Convert the list into a list of tuples, each tuple being one kinase pair/edge (since this is the format needed to 
    # add edges onto pyvis.network (network visualization package))
    edges_list2 = [tuple(x.split('.',1)) for x in edges_list]
    
    #For each treatment(column) on the file, store kinase interactions/pairs with a z-score <0 or >0, and their z-scores
    headers_in = range(1,len(headers))
    if not headers_in:
        headers_in = [1]
    if sign == 'neg':
        #Note that z-scores are converted to positive values for community detection
        edges = [[(*edges_list2[n],+ abs(networkfile[networkfile.columns[na]].to_list()[n])) for n in range(len(edges_list2)) if networkfile[networkfile.columns[na]].to_list()[n] < 0 ] for na in headers_in]
    elif sign == 'pos' :
        edges = [[(*edges_list2[n],+ networkfile[networkfile.columns[na]].to_list()[n]) for n in range(len(edges_list2)) if networkfile[networkfile.columns[na]].to_list()[n] > 0 ] for na in headers_in]
    
    #negativenetwork.add_weighted_edges_from(edges_list)
    
    #Kinases in each of the negative_edges lists
    nodes = []
    for selected_treatment in edges:
        #Keep only the kinases names (i.e. remove the z-scores) from the list of edges
        nodes_bf = [ (a,b) for a,b,c in selected_treatment]
        #Flatten list of edges tuples to list of nodes
        nodes = list(itertools.chain(*nodes_bf))
        #Remove duplicates from the list
        nodes_list = list(dict.fromkeys(nodes))
        nodes.append(nodes_list)
    
    #Convert edges and nodes list of lists into dictionary
    edges = dict(zip(treatments,edges))
    nodes = dict(zip(treatments,nodes))
    
    return treatments,edges,nodes


#1.2 - Community detection

# 1. Are all the final runs the same?
# input: edges_list is a list corresponding to a treatment in the dataset, it's contents being the kinase interactions
# in said treatment with a z_score <0 and their z_scores,in the shape of tuples, i.e. a network 
# output: A_asso is the converged association matrix, n_iter the number of iterations required for convergence and 
# labels are the converged community labels 
def converge_asso(edges_dict, sel_treatments='all'):
    
    #If no selected treatments are specified, the keys of the edges dictionary are taken as treatments
    if sel_treatments=='all':
        sel_treatments = edges_dict.keys()
    #Dictionaries with selected treatments as keys
    communities_list = dict.fromkeys(sel_treatments)
    communities_nx = dict.fromkeys(sel_treatments)
    communities_cs = dict.fromkeys(sel_treatments)

    for treatment in sel_treatments:
        #Initialize a Networkx instance
        negativenetwork = nx.Graph()
        #Add the kinase interactions & z-scores of treatment w to the network
        edges_list = edges_dict[treatment]
        negativenetwork.add_weighted_edges_from(edges_list)
        #Store the individual kinases in the network as a numpy array - to return kinases in each community at the end
        node_0 = list(negativenetwork.nodes())
        node_0_arr = np.asarray(node_0)
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
    
        #find how many communities are identified
        n_of_com = list(set(labels.flatten()))

        #List of communities of kinases as lists
        communities_list[treatment] = [list(itertools.chain(*node_0_arr[np.argwhere(labels == int(n))])) for n in n_of_com]
        #Dictionaries of communities as networkx instances
        communities = [dict.fromkeys(com,negativenetwork.subgraph(com)) for com in communities_list[treatment]]
        communities_nx[treatment] = dict(ChainMap(*communities))
        #Dictionaries of each kinase community strength
        communities_cs[treatment] = {k:v.degree([k],weight='weight')[k] for k,v in communities_nx[treatment].items()}
        
        #Add kinase community strength as attribute to nodes in the networkx dictionnaries
        for com in communities_list[treatment]:
            #Select the community strength of kinases in this community
            com_cs_dict = {kin:communities_cs[treatment][kin] for kin in com}
            #Set as node attributes in community nx instance of these kinases
            for kin in com:
                nx.set_node_attributes(communities_nx[treatment][kin], com_cs_dict, "community_strength_centrality")
    
    #Return dictionaries
    return communities_list, communities_nx, communities_cs

#1.3 - Reproducibility check

#Run converge_asso an n_checks number of times and check if the contents of the community of interest differ
def rep_check(n_checks, edges_dict, treatment, target):
    
    ##Run converge_asso n_checks times and store the contents of the identified communities each time
    all_runs_list = ['a'] * n_checks
    all_runs_nx = ['a'] * n_checks
    for n in range(n_checks):
        [all_runs_list[n],all_runs_nx[n],cs] = converge_asso(edges_dict, sel_treatments=[treatment])
  
    ##Check if there are any kinases that are not present in all the community versions
    all_runs_target_community = [ [com for com in all_runs_list[n][treatment] if target in com] for n in range(n_checks)]
    all_com = list(itertools.chain(*itertools.chain(*all_runs_target_community)))
    target_community_rem = [x for x in list(set(all_com)) if all_com.count(x)==n_checks]
    
    ##Change list of communities and dictionary on networkx instances accordingly
    
    #List
    communities_list_rem= [target_community_rem if target in com else com for com in all_runs_list[0][treatment]]
    
    #NX instances dictionary
    all_runs_nx_sel = all_runs_nx[0][treatment]
    #
    #Update remaining kinases
    community_as_nx = all_runs_nx_sel[target].subgraph(target_community_rem)
    communities_nx_rem = dict(all_runs_nx_sel, **dict.fromkeys(target_community_rem, community_as_nx))
    #
    #Assign kinases not present from all runs to X_removed
    kin_del= [kin for kin in all_runs_nx_sel[target] if kin not in target_community_rem]
    communities_nx_rem_del = dict(communities_nx_rem, **dict.fromkeys(kin_del, f"{target}<{n_checks}"))
    communities_nx_rem_del[f"{target}<{n_checks}"] = kin_del
    
    #CS dictionary
    communities_cs = {k:v.degree([k],weight='weight')[k] for k,v in communities_nx_rem.items()}
    
    #Add CS to nx instances
    for k in communities_nx_rem.keys():
        if type(communities_nx_rem_del[k]) != str:
            nx.set_node_attributes(communities_nx_rem_del[k], communities_cs, "community_strength_centrality")
        
    
    return communities_list_rem, communities_nx_rem_del, communities_cs

 #1.4 - Comparisson of CS across communities visualized in heatmap

 # input: edges_list is a list corresponding to a treatment in the dataset, it's contents being the kinase interactions
# in said treatment with a z_score <0 and their z_scores,in the shape of tuples, i.e. a network
# community is a list of the kinases in the community of interest
# output: dictionary of kinases and their community strength
def cs_comparisson(communities_nx, communities_cs, treatment_targets_sign, kin_int=None, position=0):
    #Calculate the community strength of each of the kinases in the two selected communities you are comparing
    #kin_cs = ['a'] * len(treat_int)
    #for n in range(len(treat_int)):
        #na = treat_int[n]
        # Store kinase pairs in which both kinases are in the community of interest
        #community_edges = [x for x in edges[na] if x[0] in communities[na] and x[1] in communities[na]]
        # Calculate the community strength of each of the kinases in the community and store as dictionary
        #kin_cs[n] = {kin:abs(round(sum([community_edges[n][2] for n in range(len(community_edges)) if community_edges[n][1]==kin or community_edges[n][0]==kin]), 2)) for kin in communities[na]}    

    #Specify which kinases you want to get the values of
    if kin_int == None:
        #Combine two communities lists
        kinases = set(itertools.chain(*[list(communities_nx[x[0]][x[1]] ) for x in treatment_targets_sign]))
    else:
        kinases = kin_int
    
    #Get values for each of the kinases in the combined list in each of the communities
    values = [ [communities_cs[x[0]][kin] if kin in list(communities_nx[x[0]][x[1]]) else None for kin in kinases] for x in treatment_targets_sign ]
    
    #Define heatmap axis and text
    x = list(kinases)
    #<sub></sub> for subscript, <sup></sup>
    y = [ f"{x[0]}<sup>{x[2]}</sup><sub>({x[1]})</sub><br>absolute CSC</br>" for x in treatment_targets_sign]
    z = values#[values[0],values[1]]
    
    # HEATMAP
    fig = go.Figure(data=go.Heatmap(
                    z=z,
                    x=x,
                    y=y,
                    yaxis='y2',
                    hoverongaps = False))

    fig.update_layout(xaxis=dict(domain=[0.21,1]),
                      yaxis2=dict(anchor='free', position=position, side='right'))
    
    return fig

#---------------------------------------------------------------------------------

#2 - Example 1 - Analysis and Results---------------------------------------------

###Negative edges---------------------

#Extracting networks trametinib(-), GDC0994(-), GDC0941(-), and AZD5363(-) from edges z-score dataset
#Store dataset treatments names, kinase interactions (network edges) with a z-score < 0, and their z-scores (edges weight) 
[treatments, negative_edges, negative_nodes] = dataset_info('input/PROJECT_DATASET_2.csv', rem_vs_control='yes')
#Main kinase targets of each of the treatments/inhibitors in the dataset
targets = ['MAP2K1','PIK3CA','AKT1_2','MAPK1_3']

treatments

#Networks trametinib(-) and GDC0994(-) (MEK/ERK pathway targets)
#Community detection
#Please note that the following block of code takes a few minutes to run

#Treatments in which the main target is a MEK/ERK kinase
MEK_ERK_treatments = ['Trametinib','GDC0994']
#Run community detection on these treatments negative networks
[MEK_ERK_communities, MEK_ERK_communities_nx, MEK_ERK_cs] = converge_asso(edges_dict=negative_edges, sel_treatments=MEK_ERK_treatments)

#Because we noticed that when extracting the community of interest from the trametinib(MEK)[-] network with the converge_asso
#function the contents of said community differ from run from run, for this network we run the community detection function (converge_asso)
#50 times, and keep only the nodes that are present in the community in all 50 runs of the converge_asso function
treatment = 'Trametinib'
target = 'MAP2K1'

[MEK_ERK_communities[treatment], MEK_ERK_communities_nx[treatment], MEK_ERK_cs[treatment]] = rep_check(n_checks=50,edges_dict=negative_edges,treatment=treatment,target=target)

nx.get_node_attributes(MEK_ERK_communities_nx['GDC0994']['PIK3CA'], "community_strength_centrality")['PIK3CA']

#Comparisson of community strength of communities trametinib(-)(MAP2K1) & GDC0994(-)(MAPK1_3) and visualization in heatmap
MEK_ERK_int = [('Trametinib','MAP2K1','-'),('GDC0994','MAPK1_3','-')]

fig = cs_comparisson(communities_nx=MEK_ERK_communities_nx, communities_cs=MEK_ERK_cs, treatment_targets_sign=MEK_ERK_int)

fig.show()

#Networks GDC0941(-) and AZD5363(-) (PI3K/AKT targets)
#Community detection

#Treatments in which the main target is a PIK3CA/AKT/mTOR kinase
PI3K_AKT_treatments = treatments[1:3]
#Run community detection on these treatments negative networks
[PI3K_AKT_communities, PI3K_AKT_communities_nx, PI3K_AKT_cs] = converge_asso(edges_dict=negative_edges, sel_treatments=PI3K_AKT_treatments)

#Comparisson of community strength of communities GDC0941(-)(PIK3CA) & AZD5363(-)(AKT1_2) and visualization in heatmap
PI3K_AKT_int = [('GDC0941','PIK3CA','-'),('AZD5363','AKT1_2','-')]

fig = cs_comparisson(communities_nx=PI3K_AKT_communities_nx, communities_cs=PI3K_AKT_cs, treatment_targets_sign=PI3K_AKT_int)

fig.show()

#Intersection and differences between communities trametinib(-)(MAP2K1), GDC0994(-)(MAPK1_3), GDC0941(-)(PIK3CA) and AZD5363(-)(AKT1_2)
group = com_groups_to_compare[0]

###Create communities names
com_names = [f"({x[1]}({x[0]}){x[2]})" for x in group]
    
###If comparing two communities
g = group
com1 = allcom_nx[g[0][1]][g[0][2]]
com2 = allcom_nx[g[1][1]][g[1][2]]
    
#All nodes in each community and intersection
    
com_int_and_diff = {f"{com_names[0]}_only":com1,f"{com_names[1]}_only":com2,f"{com_names[0]}{com_names[1]}":intersection(com1,com2)}
    
#Differences between communities
inter = com_int_and_diff[f"{com_names[0]}{com_names[1]}"]
com = [com1,com2]
com_int_and_diff.update({f"{com_names[n]}_only":com[n].subgraph(set(com[n].nodes()).difference(set(inter.nodes()))) for n in [0,1]})

com_int_and_diff

G = nx.path_graph(3)# or DiGraph, MultiGraph, MultiDiGraph, etc
R = nx.Graph(G)
e = nx.path_graph(2)
T = R.remove_nodes_from(e)

list(T)

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    s1 = list(chain.from_iterable(combinations(s, r) for r in range(1,len(s))))
    #s2 = ["".join(x) for x in s1]
    #s2_only = [f"{x}_only" for x in s2]
    
    #s3 = list(itertools.chain(*zip(s2,s2_only)))
    return s1

#
int_and_diff = list(itertools.chain(*[ [com_names, "".join(com_names)], [f"{com}_only" for com in com_names] ]))

int_and_diff

com_groups_to_compare = [[('Trametinib','MAP2K1'),('GDC0994','MAPK1_3')],[('GDC0941','PIK3CA'),('AZD5363','AKT1_2')],[('Trametinib','MAP2K1'),('GDC0941','PIK3CA'),('AZD5363','AKT1_2')]]
group = com_groups_to_compare[2]

#All combinations
all_int = powerset(group)
#Key
#f"({x[1]}({x[0]}){x[2]})"
#Treatment
all_int[0][0][0]
#Target
all_int[0][0][1]

#int_names = [ f"({all_int[n][0][1]}({all_int[n][0][0]}){all_int[n][0][2]})" for n in range(len(all_int))]
#all_int = {}

print(int_names)
re.sub('(', '', str(all_int[3]) )

nx.intersection_all

#Could networkx difference be useful here
#Nope, two graphs must have the same nodes
#Use intersection instead
#Unrelated, but common_neighbours could be useful for different analyses
#Finally found correct function, compose and compose_all
#Also difference and symmetric difference (nope, must have same nodes!)

#allc = [MEK_ERKallc[0],PI3K_AKTallc[0],PI3K_AKTallc[1],MEK_ERKallc[1]]
#inc = [MEK_ERKinc[0],PI3K_AKTinc[0],PI3K_AKTinc[1],MEK_ERKinc[1]]

#kinase_int = list(zip(allc,inc))

allcom_nx = dict(ChainMap(*[MEK_ERK_communities_nx, PI3K_AKT_communities_nx]))
allcom_cs = dict(ChainMap(*[MEK_ERK_cs, PI3K_AKT_cs]))

com_groups_to_compare = [[('Trametinib','MAP2K1'),('GDC0994','MAPK1_3')],[('GDC0941','PIK3CA'),('AZD5363','AKT1_2')],[('Trametinib','MAP2K1'),('GDC0941','PIK3CA'),('AZD5363','AKT1_2')]]

#R.remove_nodes_from(n for n in G if n in H)

all_int = {}
all_diff = {}
for group in com_groups_to_compare:
    
    ###Create communities names
    com_names = [f"({x[1]}({x[0]}){x[2]})" for x in group]
    
    ###Make dictionnary with all intersection and differences names
    
    
    ###If comparing two communities
    if len(group) == 2:
        g = group
        com1 = allcom_nx[g[0][1]][g[0][2]]
        com2 = allcom_nx[g[1][1]][g[1][2]]
    
        #All nodes in each community and intersection
    
        com_int_and_diff = {f"{com_names[0]}_only":com1,f"{com_names[1]}_only":com2,f"{com_names[0]}{com_names[1]}":intersection(com1,com2)}
    
        #Differences between communities
        inter = com_int_and_diff[f"{com_names[0]}{com_names[1]}"]
        com = [com1,com2]
        com_int_and_diff.update({f"{com_names[n]}_only":com[n].subgraph(set(com[n].nodes()).difference(set(inter.nodes()))) for n in [0,1]})
    
    #Intersection between all communities of interest
    #com_int = set(communities_nodes[0]).intersection(*communities_nodes)
    #Difference between each community and the intersection of all communities
    #com_diff = [set(x).difference(com_int) for x in communities_nodes]
    
    ###If comparing three communities
    #biovenn = draw_venn(list_A, list_B, list_C, subtitle="Example diagram", nrtype="abs")

    #replace = dict(zip(['x','y','z'],com_names))

    #com_int_and_diff = {multiple_rep(replace,k):v for k,v in biovenn.items()}
    
    
com_diff

#Saving the edges and nodes of each community as dataframes
#Individual communities
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

#Intersection of trametinib, GDC0941 and AZD5363
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

#-------------------------------------

###Positive edges---------------------

#Extracting networks trametinib(+), GDC0994(+), GDC0941(+), and AZD5363(+) from edges z-score dataset
#Store dataset treatments names, kinase interactions (network edges) with a z-score > 0, and their z-scores (edges weight) 
[treatments,positive_edges, positive_nodes] = dataset_info('input/PROJECT_DATASET_2.csv', sign = 'pos', rem_vs_control='yes')

treatments

#Community detection across all extracted positive edges networks
#Separate by MEK/ERK and PI3K/AKT to compare similarities in upregulated kinases
#Treatments in which the main target is a MEK/ERK kinase
MEK_ERK_treatments = ['Trametinib','GDC0994']
#Run community detection on these treatments positive networks
[MEK_ERK_pos_communities, MEK_ERK_pos_communities_nx, MEK_ERK_pos_cs] = converge_asso(edges_dict=positive_edges, sel_treatments=MEK_ERK_treatments)

MEK_ERK_pos_communities

#Treatments in which the main target is a PIK3CA/AKT/mTOR kinase
PI3K_AKT_treatments = treatments[1:3]
#Run community detection on these treatments positive networks
[PI3K_AKT_pos_communities, PI3K_AKT_pos_communities_nx, PI3K_AKT_pos_cs] = converge_asso(edges_dict=positive_edges, sel_treatments=PI3K_AKT_treatments)

PI3K_AKT_pos_communities_nx.keys()

#Comparisson of communities across network pairs

#MEK_ERK
treatments = ['Trametinib', 'GDC0994']
shared_kin_nx = nx.intersection(MEK_ERK_pos_communities_nx['Trametinib']['MAPK14'],MEK_ERK_pos_communities_nx['GDC0994']['MAPK14'] )
nodes_inf = []
for n in range(len(treatments)):
    #Kinases present in both
    shared_kin = list(shared_kin_nx)
    nodes_inf0 = list(zip(shared_kin, ['kinases present in both SCs']*len(shared_kin) ,  ['#FFDB33']* len(shared_kin) )) 
    
    #Kinases present in that network only
    not_shared_kin = [x for x in list(MEK_ERK_pos_communities_nx[treatments[n]]['MAPK14']) if x not in shared_kin]
    nodes_inf1 = list(zip( not_shared_kin, ['kinases present exclusively in this SC']*len(not_shared_kin) ,  ['#96efff']* len(not_shared_kin) )) 
    
    #Put together
    nodes_inf_sum = nodes_inf0 + nodes_inf1
    nodes_df = pd.DataFrame( nodes_inf_sum, columns=['id', 'group', 'color'])
    #add community strength
    kin_dict = {key:value for key,value in MEK_ERK_pos_cs[treatments[n]].items() if key in list(MEK_ERK_pos_communities_nx[treatments[n]]['MAPK14'])}
    kin_cs = list(zip(kin_dict.keys(), kin_dict.values() ))
    nodes_cs_df = pd.DataFrame (kin_cs, columns = ['id','community_strength'])
    nodes_df_add = pd.merge(nodes_df, nodes_cs_df, on='id')
    #save nodes as csv file
    nodes_df_add.to_csv('output_dataframes/MAPK14.{}_pos_nodes.csv'.format(treatments[n]))
    #save edges as csv file
    edges = nx.get_edge_attributes(MEK_ERK_pos_communities_nx[treatments[n]]['MAPK14'], 'weight')
    from_ed = [x[0] for x in edges.keys() ]
    to_ed = [x[1] for x in edges.keys() ]
    edges_list = list(zip(from_ed,to_ed, edges.values(),edges.values() ))
    edges_df = pd.DataFrame(edges_list, columns=['from','to','z-score','weight'])
    edges_df.to_csv('output_dataframes/MAPK14.{}_pos_edges.csv'.format(treatments[n]))

#PI3K_AKT
treatments = ['GDC0941', 'AZD5363']
shared_kin_nx = nx.intersection(PI3K_AKT_pos_communities_nx[treatments[0]]['MAPK14'],PI3K_AKT_pos_communities_nx[treatments[1]]['MAPK14'] )
nodes_inf = []
for n in range(len(treatments)):
    #Kinases present in both
    shared_kin = list(shared_kin_nx)
    nodes_inf0 = list(zip(shared_kin, ['kinases present in both SCs']*len(shared_kin) ,  ['#FFDB33']* len(shared_kin) )) 
    
    #Kinases present in that network only
    not_shared_kin = [x for x in list(PI3K_AKT_pos_communities_nx[treatments[n]]['MAPK14']) if x not in shared_kin]
    nodes_inf1 = list(zip( not_shared_kin, ['kinases present exclusively in this SC']*len(not_shared_kin) ,  ['#96efff']* len(not_shared_kin) )) 
    
    #Put together
    nodes_inf_sum = nodes_inf0 + nodes_inf1
    nodes_df = pd.DataFrame( nodes_inf_sum, columns=['id', 'group', 'color'])
    #add community strength
    kin_dict = {key:value for key,value in PI3K_AKT_pos_cs[treatments[n]].items() if key in list(PI3K_AKT_pos_communities_nx[treatments[n]]['MAPK14'])}
    kin_cs = list(zip(kin_dict.keys(), kin_dict.values() ))
    nodes_cs_df = pd.DataFrame (kin_cs, columns = ['id','community_strength'])
    nodes_df_add = pd.merge(nodes_df, nodes_cs_df, on='id')
    #save nodes as csv file
    nodes_df_add.to_csv('output_dataframes/MAPK14.{}_pos_nodes.csv'.format(treatments[n]))
    #save edges as csv file
    edges = nx.get_edge_attributes(PI3K_AKT_pos_communities_nx[treatments[n]]['MAPK14'], 'weight')
    from_ed = [x[0] for x in edges.keys() ]
    to_ed = [x[1] for x in edges.keys() ]
    edges_list = list(zip(from_ed,to_ed, edges.values(),edges.values() ))
    edges_df = pd.DataFrame(edges_list, columns=['from','to','z-score','weight'])
    edges_df.to_csv('output_dataframes/MAPK14.{}_pos_edges.csv'.format(treatments[n]))

#-------------------------------------

#---------------------------------------------------------------------------------

#3 - Example 2 - Analysis and Results---------------------------------------------

###Negative edges---------------------

#Store dataset treatments names, kinase interactions (network edges) with a z-score < 0, and their z-scores (edges weight) 
[treatments2,negative_edges2,negative_nodes2] = dataset_info('input/example2_edges_zscores.csv')
#Main kinase targets of each of the treatments/inhibitors in the dataset
targets2 = ['MTOR']

treatment = ['Rapamycin_vs_Control']
#Store contents of communities identified in each network
[rap_ER_communities, rap_ER_communities_nx, rap_ER_cs] = converge_asso(edges_dict=negative_edges2, sel_treatments=treatment)
#Reproducibiliyyt check
treatment = 'Rapamycin_vs_Control'
target = 'MTOR'
#[rap_ER_communities[treatment], rap_ER_communities_nx[treatment], rap_ER_cs[treatment]] = rep_check(n_checks=50,edges_dict=negative_edges2,treatment=treatment,target=target)

treatments = list(PI3K_AKT_pos_communities_nx.keys())
PI3K_AKT_shared = nx.intersection(PI3K_AKT_communities_nx[treatments[0]]['PIK3CA'],PI3K_AKT_communities_nx[treatments[1]]['AKT1_2'] )
shared_kin_nx = nx.intersection(rap_ER_communities_nx['Rapamycin_vs_Control']['MTOR'], PI3K_AKT_shared)
#Also shared with either MAP2K1 or MAPK1_3
MEK_ERK = nx.intersection(MEK_ERK_communities_nx['Trametinib']['MAP2K1'], MEK_ERK_communities_nx['GDC0994']['MAPK1_3'])
MEK_ERK_shared = nx.intersection(rap_ER_communities_nx['Rapamycin_vs_Control']['MTOR'], MEK_ERK)

# Save name of each of the communities of interest to show in the graph visualizations (see below)

PI3K_AKT_int = [('Rapamycin_vs_Control','MTOR','-'), ('GDC0941','PIK3CA','-'),('AZD5363','AKT1_2','-')]
sup_st = '$^{'
sup_end = '}$'
sub_st = '$_{'
sub_end = '}$'
PI3K_AKT_names = [ f"{x[0]}{sup_st}{x[2]}{sup_end}{sub_st}({x[1]}){sub_end}" for x in PI3K_AKT_int]

MEK_ERK_int = [('Rapamycin_vs_Control','MTOR','-'), ('trametinib','MAP2K1','-'),('GDC0994','MAPK1_3','-')]
MEK_ERK_names = [ f"{x[0]}{sup_st}{x[2]}{sup_end}{sub_st}({x[1]}){sub_end}" for x in MEK_ERK_int]

##Graph visualizations

treatments = ['Rapamycin_vs_Control']
com = 'MTOR'
treat = 'MCF7_AML_com_comparisson'
ext='pdf'
for n in range(len(treatments)):
    net = rap_ER_communities_nx[treatments[n]][com]
    #pigraphviz instance from networkx graph
    A = nx.nx_agraph.to_agraph(net)
    #Map nodes weights to colors
    color_lookup = nx.get_node_attributes(net, 'community_strength_centrality')
    min_value, *_, max_value = sorted(color_lookup.values())
    norm = mpl.colors.Normalize(vmin=min_value, vmax=max_value, clip=True)
    mapper = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.GnBu) #Blues)
    other_mapper = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Oranges)#YlOrBr)
    #Map edges weights to darker or lighter grey based on weight
    ed_color_lookup = nx.get_edge_attributes(net, 'weight')
    min_value, *_, max_value = sorted(ed_color_lookup.values())
    norm = mpl.colors.Normalize(vmin=min_value, vmax=max_value, clip=True)
    edges_mapper = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Greys) #Blues)
    for edge in ed_color_lookup.keys():
        ed = A.get_edge(edge[0],edge[1])
        ed.attr['color'] = mcolors.to_hex(edges_mapper.to_rgba(ed_color_lookup[edge]))
        ed.attr['penwidth'] = ed_color_lookup[edge]*1
    #Draw nodes above edges
    A.graph_attr["strict"] = False
    A.graph_attr["overlap"] = False
    A.graph_attr["splines"] = 'compound'
    #Color by group
    A.node_attr['style']='filled'
    ##Nodes present in this community but not the other(s)
    #unique_nodes = [node for node in list(net) if node not in shared_kin_nx]
    for node in list(MEK_ERK_shared):
        n = A.get_node(node)
        n.attr['fontsize'] = 22
        n.attr['fontname'] = "helvetica bold"
        #n.attr['fillcolor']="#9bd9f2"
        rgba = other_mapper.to_rgba(color_lookup[node])
        n.attr['fillcolor']= mcolors.to_hex(rgba)
        #Change font color to white if background is too dark
        if (rgba[0]*0.299 + rgba[1]*0.587 + rgba[2]*0.114) < font_threshold:
            n.attr['fontcolor']= '#ffffff'
    ##Nodes present in both communities being compared
    for node in list(shared_kin_nx):
        n = A.get_node(node)
        n.attr['fontsize'] = 22
        n.attr['fontname'] = "helvetica bold"
        #n.attr['fillcolor']="#9bd9f2"
        rgba = mapper.to_rgba(color_lookup[node])
        n.attr['fillcolor']= mcolors.to_hex(rgba)
        #Change font color to white if background is too dark
        if (rgba[0]*0.299 + rgba[1]*0.587 + rgba[2]*0.114) < font_threshold:
            n.attr['fontcolor']= '#ffffff'
    for node in list(net):
        if node not in list(MEK_ERK_shared) and node not in list(shared_kin_nx):
            n = A.get_node(node)
            n.attr['fontsize'] = 22
            n.attr['fontname'] = "helvetica bold"
            
    #Highlight MEK/ERK
    for node in ['ARAF','MAP2K1','MAPK1_3']:
        if node in list(net):
            n = A.get_node(node)
            n.attr['shape']='square'
            n.attr['width']=1.5
            n.attr['fixedsize']='true' 
            #n.attr['penwidth']= 3
            #n.attr['color']="#FFF300"
            #n.attr['color']="#4AF1F2"
            #n.attr['fillcolor']="#9bd9f2"
            #n.attr['fillcolor']= mcolors.to_hex(mapper.to_rgba(color_lookup[node]))
    #Highlight PIK3CA/AKT
    for node in ['PIK3CA','AKT1_2','MTOR']:
        if node in list(net):
            n = A.get_node(node)
            n.attr['shape']='diamond'
            n.attr['width']=1.5
            n.attr['height']=2
            n.attr['fixedsize']='true'
            #n.attr['penwidth']= 3
            #n.attr['color']="#FF6EFF"
            #n.attr['fillcolor']="#9bd9f2"
            #n.attr['fillcolor']= mcolors.to_hex(mapper.to_rgba(color_lookup[node]))
    
    
    ##Multi-panel colorbar
        
    #Add colorbar legend for unique nodes 
    fig, ax = plt.subplots(1, 1)
    fraction = 1
    plt.colorbar(other_mapper, ax=ax, pad=.05, extend='both', fraction=fraction)
    ax.axis('off')
    ax.set_title(f"absolute CSC of kinases \n present across \n {MEK_ERK_names[0]} \n {MEK_ERK_names[1]} and {MEK_ERK_names[2]}", x=-0.3)
    name = f"output_visualizations/{treat}_pos_networkviz_static_unique_nodes_colorbar.{ext}"
    plt.savefig(name, bbox_inches="tight",format=f"{ext}")
    plt.show()
    
    #Add colorbar legend for shared nodes
    fig, ax = plt.subplots(1, 1)
    fraction = 1
    plt.colorbar(mapper, ax=ax, pad=.05, extend='both', fraction=fraction)
    ax.axis('off')
    ax.set_title(f"absolute CSC of kinases \n present across \n {PI3K_AKT_names[0]} \n {PI3K_AKT_names[1]} and  {PI3K_AKT_names[2]}", x=-0.3)
    name = f"output_visualizations/{treat}_pos_networkviz_static_shared_nodes_colorbar.{ext}"
    plt.savefig(name, bbox_inches="tight",format=f"{ext}")
    plt.show()
    
    #Add shapes legend
    fig, ax = plt.subplots(1, 1)
    name = f"output_visualizations/MAPK14.{treat}_pos_networkviz_static_shared_nodes_colorbar.{ext}"
    a = ax.scatter([0], [0], marker='D', facecolors='none', edgecolors='black', label= 'PI3K/AKT/mTOR')
    b = ax.scatter([0], [0], marker='s', facecolors='none', edgecolors='black', label= 'ARAF/MAP2K1/MAPK1_3' )
    # adjust the figure size as necessary
    fig_leg = plt.figure(figsize=(1,1))
    ax_leg = fig_leg.add_subplot(111)
    # add the legend from the previous axes
    ax_leg.legend(*ax.get_legend_handles_labels(), loc='center', title='canonical pathways')
    # hide the axes frame and the x/y labels
    ax_leg.axis('off')
    ax.set_title(title_shapes, x=-0.3)
    name = f"output_visualizations/{treat}_pos_networkviz_static_pathways_legend.{ext}"
    fig_leg.savefig(name, bbox_inches="tight",format=f"{ext}")
    fig_leg.show()

    #Add colorbar legend for edges
    fig, ax = plt.subplots(1, 1)
    fraction = 1
    plt.colorbar(edges_mapper, ax=ax, pad=.05, extend='both', fraction=fraction)
    ax.axis('off')
    ax.set_title('absolute edges z-score \n (KSEA)', x=-0.3)
    name = f"output_visualizations/{treat}_pos_networkviz_static_edges_colorbar.{ext}"
    plt.savefig(name, bbox_inches="tight",format=f"{ext}")
    plt.show()
    
    
    #Draw graph
    A.layout()
    title = f"output_visualizations/{treat}_pos_networkviz_static"
    A.draw(f"{title}.{ext}", prog='sfdp',format=f"{ext}")
    #os.system(f"pdftoppm -png -r 300  {title}.{ext} {title}")

## Heatmap of kinases CS across communities:
# rapamycin(-)(mTOR)
# AZD5363(-)(AKT)
# GDC0941(-)(PI3K)

import plotly.io as pio
import kaleido
## Comparative heatmap

#Combine dictionnaries
all_com = dict(ChainMap(*[rap_ER_communities_nx, PI3K_AKT_communities_nx]))
all_com_cs = dict(ChainMap(*[rap_ER_cs, PI3K_AKT_cs]))
for net in all_com.keys():
    all_com[net] = {kin:all_com[net][kin] for kin in all_com[net].keys() if kin in shared_kin_nx}
    all_com_cs[net] = {kin:all_com_cs[net][kin] for kin in all_com_cs[net].keys() if kin in shared_kin_nx}

#Figure
all_com['Rapamycin'] = all_com.pop('Rapamycin_vs_Control')
all_com_cs['Rapamycin'] = all_com_cs.pop('Rapamycin_vs_Control')
com_int = [('Rapamycin','MTOR','-'), ('GDC0941','PIK3CA','-'),('AZD5363','AKT1_2','-')]

fig = cs_comparisson(communities_nx=all_com, communities_cs=all_com_cs, treatment_targets_sign=com_int, kin_int=shared_kin_nx, position=0)
#fig1 = fig.gcf()
#fig.write_image(file='output_visualizations/MCF7_AML_com_comparisson.pdf', format='pdf', engine='kaleido')
fig.update_layout(
    font=dict(
        size=16,
    )
)
fig.show()

#-------------------------------------

###Positive edges---------------------

#Store dataset treatments names, kinase interactions (network edges) with a z-score < 0, and their z-scores (edges weight) 
[treatments2,positive_edges2, positive_nodes2] = dataset_info('input/example2_edges_zscores.csv', sign='pos')

treatment = ['Rapamycin_vs_Control']
#Store contents of communities identified in each network
[rap_ER_pos_communities, rap_ER_pos_communities_nx, rap_ER_pos_cs] = converge_asso(edges_dict=positive_edges2, sel_treatments=treatment)
#Reproducibiliyyt check
treatment = 'Rapamycin_vs_Control'

rap_ER_pos_communities['Rapamycin_vs_Control']

#-------------------------------------

#---------------------------------------------------------------------------------





