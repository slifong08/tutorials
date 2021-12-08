#!/usr/bin/env python
# coding: utf-8

# In[1]:


import configparser
import os
import subprocess
import sys


# # load config

# In[2]:


DEV = True

if DEV is True:
    BASE_PATH = "/dors/capra_lab/users/fongsl/tutorials/using_config_ini/"
    configfile_name = os.path.join(BASE_PATH, "bin", "config_gwas.ini") # name the file
else:
    configfile_name = sys.argv[1]
    
config = configparser.ConfigParser()
config.read(configfile_name)


# # load parameters 

# In[3]:


GWAS_CAT_CLEAN_BED = config["FILE"]["GWAS_CAT_CLEAN_BED"]  # read
GWAS_CAT_CLEAN_LD = config["FILE"]["GWAS_CLEAN_LD"]  # write
ASSOC2LD = config["BIN"]["LD_ASSOC_SCRIPT"]  # bin

POP = config["PARAMS"]["POP"]  # param


# # Functions 

# In[4]:


def run_assoc2SNP(assoc2ldsnp, pop, gwas_cat_clean_bed, gwas_cat_clean_expanded):
    
    cmd = f"python {assoc2ldsnp} {pop} {gwas_cat_clean_bed} > {gwas_cat_clean_expanded}"

    if os.path.exists(gwas_cat_clean_expanded) is False:
        print("LD expanding")
        subprocess.call(cmd, shell = True) 
        
    else:
        print("already expanded")


# # Main 

# In[5]:


def main(argv):
    
    run_assoc2SNP(ASSOC2LD, POP, GWAS_CAT_CLEAN_BED, GWAS_CAT_CLEAN_LD)
    
if __name__ == "__main__":
    main(sys.argv[1:])


# # write notebook as .py

# In[6]:


NAME= "LD_expand_GWAS"
cmd = f"jupyter nbconvert --to python {NAME}.ipynb"
subprocess.call(cmd, shell = True)

