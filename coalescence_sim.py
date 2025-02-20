#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#pip install msprime -q

import msprime
import tskit

import numpy as np
import random as rd
import math

import pandas as pd
import pickle

import os


# 
# # Methods: Simulate unlinked loci, compute allele frequencies, and generate genotype matrices for a split or migration model

# In[ ]:


### Simulate an independent coalescent tree for an N-deme split model
def sim_split_singtree(Ne, n_samp, d, t, subtype):

  n = Ne            # int: local population size
  n_samp = n_samp   # int: number of individuals sampled per deme (number of chromosomes per deme = 2*n_samp)
  N = d             # int: number of demes
  t = t             # float: split time (overlying for the balanced split model)
  n_anc = Ne        # int: number of individuals in the ancestral population (here, n_anc=Ne for basic models)

  if (N == 1):
    demography = None
    dicts = n_samp
  else:
    demography = msprime.Demography()
    for i in range(N):
      demography.add_population(initial_size=n)

    names = [f"p{i}" for i in range(N)]
    for i in range(N):
      demography[i].name = names[i]
    dicts = {}
    keys = names
    for i in keys:
      dicts[i] = n_samp

    if subtype == "star":
      demography.add_population(initial_size=n_anc)
      demography[N].name = f"p{N}"
      demography.add_population_split(time=t, derived=names, ancestral=f"p{N}")
    elif subtype == "balanced":
      if ((N & (N-1) == 0) and N != 0):
        for i in range(N, 2*N-1):
          demography.add_population(initial_size=n_anc)
          demography[i].name = f"p{i}"
        
        # initialize a split time list
        t_list = np.zeros((N-1,), dtype=float)
        t_0 = t
        
        def binary_split(arr, left, right, value):
          if right - left == 1:
            arr[left] = value
            return
          mid = math.ceil((left + right) / 2)
          arr[left:mid] = value
          binary_split(arr, mid, right, value+t_0)
        
        binary_split(t_list, 0, N-1, t)

        # add sequential split events
        for j in range(N, 2*N-1):
          demography.add_population_split(time=t_list[j-N], derived=[f"p{2*(j-N)}", f"p{2*(j-N)+1}"], ancestral=f"p{j}")
      else:
        print("Error: a balanced split tree can only be generated when the number of demes is a power of 2.")
    elif subtype == "caterpillar":
      for i in range(N, 2*N-1):
        demography.add_population(initial_size=n_anc)
        demography[i].name = f"p{i}"
        if i == N:
          demography.add_population_split(time=(i-N+1)*t, derived=["p0", "p1"], ancestral=f"p{i}")
        else:
          demography.add_population_split(time=(i-N+1)*t, derived=[f"p{(i-1)}", f"p{i-N+1}"], ancestral=f"p{i}")


  ts = msprime.sim_ancestry(
      samples=dicts,
      demography=demography,
      recombination_rate=0,     # totally linked, only produce one tree
      sequence_length=1*10**6,
      random_seed = None)

  return ts


# In[ ]:


### Simulate an independent coalescent tree for an N-deme migration model

def sim_migration_singtree(Ne, n_samp, d, m, subtype):

  n = Ne            # int: local population size
  n_samp = n_samp   # int: number of individuals sampled per deme (number of chromosomes per deme = 2*n_samp)
  N = d             # int: number of demes
  m = m             # float: overall migration rate in an island model,
                    # or migration rate going out from one deme to either of the two adjacent demes in a stepping-stone model

  if N == 1:
    demography = None
    dicts = n_samp
  else:
    if subtype == "island":
      demography = msprime.Demography.island_model(initial_size=[n]*N, migration_rate=m/(N-1))
    elif subtype == "circular": # circular stepping-stone
      demography = msprime.Demography.stepping_stone_model(initial_size=[n]*N, migration_rate=m/2, boundaries=False)

    names = [f"p{i}" for i in range(N)]
    for i in range(N):
      demography[i].name = names[i]
    dicts = {}
    keys = names
    for i in keys:
      dicts[i] = n_samp


  ts = msprime.sim_ancestry(
          samples=dicts,
          demography=demography,
          recombination_rate=0,     # totally linked, only produce one tree
          sequence_length=1*10**6,
          random_seed=None)


  return ts


# In[ ]:


#helper function: recursively calculate the number of terminal nodes in each population that can be reached from an internal node
#u is the node of interest. Computation for the root node returns computation for all the nodes

def bl_spectrum(tree):
  counts = dict()
  rev_counts = dict()

  u = tree.root
  n = d

  def terminalnodecounts(u):
    attr = {i: 0 for i in range(0, n)}

    if not tree.children(u):
      attr[tree.population(u)] = 1
    else:
      attr_left = terminalnodecounts(tree.left_child(u))
      attr_right = terminalnodecounts(tree.right_child(u))
      list1 = list(attr_left.keys())
      list2 = list(attr_right.keys())
      shared_keys = set(list1).intersection(list2)
      if not shared_keys:
        attr = {**attr_left, **attr_right}
      else:
        new = dict()
        for i in shared_keys:
          value = attr_left[i] + attr_right[i]
          new[i] = value
        attr = {**attr_left, **attr_right, **new}

    counts[u] = tuple(attr.values())

    return attr

  _ = terminalnodecounts(u)

  for k, v in counts.items():
    rev_counts[v] = rev_counts.get(v, []) + [k]

  return rev_counts


# In[ ]:


### Store the observed branch lengths in the spectrum into a B-set (as dictionary)

def getBset(tree, rev_counts):
  res = dict()

  for k, v in rev_counts.items():
    sum_bl = 0
    for i in range(0, len(v)):
      sum_bl = sum_bl + tree.branch_length(v[i])
    res[k] = sum_bl

  res["total branch length"]= tree.total_branch_length

  return res


# In[ ]:


### Calculate demographic parameters for models with a specific FST
### split time (t) for split models, migraiton rate (m) for migration models

def getParam(subtype, d, Ne, FST):
  formulas = {
        "star": FST/(1-FST) * 2*Ne*d/(d-1),
        "balanced": FST/(1-FST) * 2*Ne*d/(d*math.log2(d)-d+1),
        "caterpillar": FST/(1-FST) * 6*Ne*d/((d-1)*(2*d-1)),
        "island": (1-FST)/FST * ((d-1)**2)/(4*Ne*d**2),
        "circular": (1-FST)/FST * ((d**2)-1)/(24*Ne*d),
  }

  return formulas.get(subtype, "Invalid parameter")


# In[ ]:


### Run multiple simulations and generate a meta branch length spectrum of independent loci

def bl_spectrum_indloci(n_sim):
  meta_res = dict()

  for i in range(n_sim):
    if i > 0 and i % 100 == 0:
      #temp = open('/content/drive/MyDrive/temp_0420/temp'+str(i)+'_meta_res.txt', 'wb')
      #pickle.dump(meta_res, temp)
      #temp.close()
      print(str(i)+ "saved")

    if event == "split":
      t = getParam(subtype, d, Ne, FST)
      ts = sim_split_singtree(Ne, n_samp, d, t, subtype)
    elif event == "migration":
      m = getParam(subtype, d, Ne, FST)
      ts = sim_migration_singtree(Ne, n_samp, d, m, subtype)

    tree = ts.first()
    rev_counts = bl_spectrum(tree)
    res = getBset(tree, rev_counts)

    list1 = list(meta_res.keys())
    list2 = list(res.keys())
    shared_keys = set(list1).intersection(list2)
    if not shared_keys:
      meta_res = {**meta_res, **res}
    else:
      new = dict()
      for i in shared_keys:
        value = meta_res[i] + res[i]
        new[i] = value
      meta_res = {**meta_res, **res, **new}

  return meta_res


# In[ ]:


### Generate SFS dictionary/matrix and/or marginal SFS matrix from the meta branch length spectrum

def sfs(B, add_matrix: bool = False, add_marginal: bool = False):
  sfs = B.copy()
  total = sfs['total branch length']
  sfs = {k: v / total for k, v in sfs.items()}
  del sfs['total branch length']
  res = [sfs]

  n_tip = 2*n_samp

  # create an SFS matrix if needed （for 2-Deme models only)
  if add_matrix == True:
    sfs_mtx = np.zeros((n_tip+1, n_tip+1))
    for k, v in sfs.items():
      i = k[0]
      j = k[1]
      sfs_mtx[i,j] = v
    res.append(sfs_mtx)

  # create a marginal SFS matrix if needed
  if add_marginal == True:
    msfs_mtx = np.zeros((n_tip, d))
    for deme in range(d):
      for k,v in sfs.items():
        if k[deme]>0 :
          msfs_mtx[k[deme]-1,deme] += v
    res.append(msfs_mtx)

  return res


# In[ ]:


### Construct the meta-population genotype matrix according to SFS

def genMtx(sfs, n_allele):
  p = sfs
  n_tip = 2*n_samp

  # generate allele frequency vector
  freq_vec = np.random.multinomial(n_allele, list(p.values()), size=1)
  p_freq = dict()
  i = 0
  for k in p.keys():
    p_freq[k] = freq_vec[0][i]
    if freq_vec[0][i] == 0:
      del p_freq[k]
    i = i+1
  #print(p_freq)

  # use a list of d (default: 10) n_allele by n_tip (default: 200000x200) array to store matrix per population
  res = [np.zeros((sum(p_freq.values()), n_tip)) for _ in range(d)]
  cut = 0
  for k, v in p_freq.items():
    pre = cut
    cut += v
    for x in range(d):
      n_allele1 = k[x]
      tmp = np.zeros((v, n_tip))
      for n in range(v):
        a = rd.sample(range(n_tip), n_allele1)
        b = np.zeros(n_tip)
        for i in a: b[i] = 1
        tmp[n, :] = b
      res[x][pre:cut, :] = tmp
  #print(res)

  # get a diploid meta-population genotype matrix
  res = [x[:, :n_samp] + x[:, n_samp:2*n_samp] for x in res]
  geno_mtx = np.hstack(res)

  return geno_mtx


# # Simulation starts here

# In[ ]:


### Initialization
event = "split"
subtype = "star"
d = 2
FST = 0.1


# In[ ]:


### Simulation
Ne = 1000
n_samp = int(2000)
n_allele = int(20000)

B = bl_spectrum_indloci(n_sim=5000)
sfs_res = sfs(B, add_matrix=False, add_marginal=True)
geno_mtx = genMtx(sfs_res[0], n_allele=n_allele)

### Saving
folder_name = f"fst{FST}"
base_path = '/Users/BidYD/Desktop/F/fst/input/'
folder_path = os.path.join(base_path, folder_name)

if not os.path.exists(folder_path):
    os.makedirs(folder_path)

file_name = f"geno_{d}D_{event}_{subtype}_test.txt"
file_path = os.path.join(folder_path, file_name)

np.savetxt(file_path, geno_mtx, fmt='%i')


# In[ ]:




