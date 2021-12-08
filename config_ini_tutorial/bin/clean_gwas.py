#!/usr/bin/env python
# coding: utf-8

# In[1]:


import configparser
import os
import pandas as pd
import sys
import subprocess


# #  load config

# In[8]:


DEV = True

if DEV is True:
    BASE_PATH=os.getcwd()
    configfile_name = os.path.join(BASE_PATH, "config_gwas.ini") # name the file

else:
    configfile_name = sys.argv[1]
    
config = configparser.ConfigParser()
config.read(configfile_name)


# #  load parameters

# In[3]:


GWAS_CAT = config["FILE"]["GWAS_CAT"]  # read

GWAS_CAT_CLEAN = config["FILE"]["GWAS_CAT_CLEAN"]  # write
GWAS_CAT_CLEAN_BED = config["FILE"]["GWAS_CAT_CLEAN_BED"]  # write

SIG_PVAL = config.getfloat("PARAMS", "SIG_PVAL")


# # Functions

# In[4]:


def clean_gwas_catalog(gwas_cat_file, sig_pval, gwas_cat_clean_file, gwas_cat_clean_bed):
    
    # open the file
    df = pd.read_csv(gwas_cat_file, sep = '\t', nrows = 500)  

    # filter for significant P-values
    clean = df.loc[df["P-VALUE"].astype(float) < sig_pval]  

    # write filtered table
    clean.to_csv(gwas_cat_clean_file, sep = '\t', index = False)  
    
    print("cleaned project file", df.shape[0], "to", clean.shape[0])
    
    # make a .bed-like dataframe
    clean_coor = clean[["#CHR", "START", "END", "SNPS"]].drop_duplicates()
    
    # write bed coordinates
    clean_coor.to_csv(gwas_cat_clean_bed, sep ='\t', index = False)


# # Main

# In[5]:


def main(argv):
    if os.path.exists(GWAS_CAT_CLEAN_BED) is False:

        clean_gwas_catalog(GWAS_CAT, SIG_PVAL, GWAS_CAT_CLEAN, GWAS_CAT_CLEAN_BED)

    else:
        print("already cleaned this gwas file")

if __name__ == "__main__":
    main(sys.argv[1:])


# # write notebook as .py

# In[6]:


NAME= "clean_gwas"
cmd = f"jupyter nbconvert --to python {NAME}.ipynb"
subprocess.call(cmd, shell = True)


# In[ ]:




