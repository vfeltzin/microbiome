import matplotlib.pyplot as plt
import numpy as np

def hist_otu_read_abundance(otu_dict):
    otu_list = []
    for item in otu_dict.values():
        otu_list.append(len(item))
    otuhist = plt.hist(otu_list,bins=[0,2,10,100,1000,10000,100000])
    plt.gca().set_xscale("log")
    plt.show()
    return(np.histogram(otu_list,bins=[0,2,10,100,1000,10000,100000]))

def plot_otu_read_abundance(otu_dict):
    otu_list = []
    for item in otu_dict.values():
        otu_list.append(len(item))
    otuhist = np.histogram(otu_list,bins=[1,10,100,1000,10000,100000])
    plt.scatter(otuhist[0],[1,10,100,1000,10000,100000])
    plt.gca().set_xscale("log")
    plt.show()
    return(otuhist)