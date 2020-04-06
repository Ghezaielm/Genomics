#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np 
import matplotlib.pyplot as plt 
from Bio import Entrez
import yaml
import pandas as pd 
from sklearn.decomposition import PCA


# In[2]:


class covid_genome_analysis(): 
    def __init__(self): 
        self.yaml = 'metadata.yaml'
        self.accessions = []
        self.dates = []
        self.countries = []
        self.metadata = pd.DataFrame()
        self.genomes = {}
    def get_metadata(self): 
        # Fetch only complete genomes
        with open(self.yaml) as f:
            docs = yaml.load_all(f, Loader=yaml.FullLoader)
            for doc in docs:
                for category, content in doc.items():
                    if category == "genbank-sequences": 
                        for genome in content:
                            if genome["gene-region"]=="complete": 
                                self.accessions+=[genome["accession"]]
                                self.dates+=[genome["collection-date"]]
                                self.countries+=[genome["locality"]["country"]]

        self.metadata["date"]=self.dates
        self.metadata["country"]=self.countries
        self.metadata["accession"]=self.accessions 
        print(self.metadata.head())
    
    def get_genomes(self): 
        Entrez.email = "Stay@Home.now"
        self.genomes = {accession:[] for accession in self.metadata['accession']}
        dic = {"A":1,"T":2,"G":3,"C":4}
        for idx, accession in enumerate(self.genomes):
            query = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta",
                                 seq_start=21000, seq_stop=26000)
            for line in query: 
                if ">" not in line: 
                    self.genomes[accession]+=[dic[nt.upper()] 
                                              for nt in line 
                                              if nt.upper() in dic]
        ref = self.genomes["MN908947"]
        for idx,genome in enumerate(self.metadata["accession"][1:]):
            next_g = self.genomes[genome] 
            print(genome)
            if len(next_g)==len(ref):
                shift = [i for i in range(1,300)]
                shift = [-i for i in shift][::-1]+[0]+shift

                for i in zip(shift,[np.corrcoef(ref,np.roll(next_g,shift=i))[0][1] for i in shift]):
                    if i[1]>0.8:
                        self.genomes[genome]=np.roll(next_g,shift=i[0])[563:-616]
                        print(idx)
            else:
                self.genomes.pop(genome)
        self.metadata = self.metadata[self.metadata["accession"].isin(self.genomes)]
o = covid_genome_analysis()
o.get_metadata()
o.get_genomes()
o.align_spike_s()


# In[22]:


class data_processing(covid_genome_analysis): 
    
    def __init__(self):
        self.genomes = o.genomes
        self.metadata = o.metadata
        
    def heat_map(self):
        import seaborn as sns
        sns.set()
        data = pd.DataFrame([self.genomes[i] for i in self.genomes]).dropna(axis=1)
        data.index = [i for i in self.genomes]

    def get_GC_rate(self):
        GC = []
        for genome in self.genomes:
            curr = self.genomes[genome]
            A_count = [i for i in curr if i==1]
            T_count = [i for i in curr if i==2]
            G_count = [i for i in curr if i==3]
            C_count = [i for i in curr if i==4]
            GC_rate = (len(G_count)+len(C_count))/len(curr)
            print(GC_rate)
            GC.append(GC_rate)
        f = plt.figure(figsize=(25,10))
        ax= f.add_subplot(111)
        ax.scatter(self.metadata["country"],GC)
        ax.tick_params(axis='x', rotation=45)
        ax.grid(b=None)
        
u = data_processing()
u.heat_map()
u.get_GC_rate()


# In[44]:




