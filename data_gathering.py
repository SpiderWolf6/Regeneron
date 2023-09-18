#Import necessary packages
import io
import re
import time
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from bs4 import BeautifulSoup
import numpy as np
import json
from urllib.request import urlopen
from selenium import webdriver
from selenium.webdriver import ActionChains
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from flask import Flask, render_template, request, redirect, session, jsonify
from flask_cors import CORS
from flask_restful import Api, Resource, reqparse
from waitress import serve

#Create three empty dictionaries for cell membrance, cell surface, and ECM proteins
mem_dict = {"Protein Name": [], "UniPROT ID": [], "Signal Sequence": [], "Expressions Scores (Highest)": [],
              "Expressions Scores (2nd Highest)": [], "Expressions Scores (3rd Highest)": [],
              "Expressions Scores (Lowest)": [], "Expressions Scores (2nd Lowest)": [],
              "Expressions Scores (3rd Lowest)": []}

surf_dict = {"Protein Name": [], "UniPROT ID": [], "Signal Sequence": [], "Expressions Scores (Highest)": [],
              "Expressions Scores (2nd Highest)": [], "Expressions Scores (3rd Highest)": [],
              "Expressions Scores (Lowest)": [], "Expressions Scores (2nd Lowest)": [],
              "Expressions Scores (3rd Lowest)": []}

ecm_dict = {"Protein Name": [], "UniPROT ID": [], "Signal Sequence": [], "Expressions Scores (Highest)": [],
              "Expressions Scores (2nd Highest)": [], "Expressions Scores (3rd Highest)": [],
              "Expressions Scores (Lowest)": [], "Expressions Scores (2nd Lowest)": [],
              "Expressions Scores (3rd Lowest)": []}


#Filtering, reorganizing, and splicing each database using the function below
def filter_sort(temp_df):
    collection1 = []
    for ind in temp_df.index:
        var = temp_df["Expressions Scores (Highest)"][ind]
        if var.find("See") != -1:
            x = re.findall("\d+\.\d+", var)
            start = var.find(x[0])
            y = var[:start]
            collection1.append(y)
            temp_df["Expressions Scores (Highest)"][ind] = x[0]
        else:
            collection1.append("null")
    temp_df["Expressions Scores (Highest) - Tissue Type"] = collection1

    collection2 = []
    for ind in temp_df.index:
        var = temp_df["Expressions Scores (2nd Highest)"][ind]
        if var.find("See") != -1:
            x = re.findall("\d+\.\d+", var)
            start = var.find(x[0])
            y = var[:start]
            collection2.append(y)
            temp_df["Expressions Scores (2nd Highest)"][ind] = x[0]
        else:
            collection2.append("null")
    temp_df["Expressions Scores (2nd Highest) - Tissue Type"] = collection2

    collection3 = []
    for ind in temp_df.index:
        var = temp_df["Expressions Scores (3rd Highest)"][ind]
        if var.find("See") != -1:
            x = re.findall("\d+\.\d+", var)
            start = var.find(x[0])
            y = var[:start]
            collection3.append(y)
            temp_df["Expressions Scores (3rd Highest)"][ind] = x[0]
        else:
            collection3.append("null")
    temp_df["Expressions Scores (3rd Highest) - Tissue Type"] = collection3

    collection4 = []
    for ind in temp_df.index:
        var = temp_df["Expressions Scores (Lowest)"][ind]
        if var.find("See") != -1:
            x = re.findall("\d+\.\d+", var)
            start = var.find(x[0])
            y = var[:start]
            collection4.append(y)
            temp_df["Expressions Scores (Lowest)"][ind] = x[0]
        else:
            collection4.append("null")
    temp_df["Expressions Scores (Lowest) - Tissue Type"] = collection4

    collection5 = []
    for ind in temp_df.index:
        var = temp_df["Expressions Scores (2nd Lowest)"][ind]
        if var.find("See") != -1:
            x = re.findall("\d+\.\d+", var)
            start = var.find(x[0])
            y = var[:start]
            collection5.append(y)
            temp_df["Expressions Scores (2nd Lowest)"][ind] = x[0]
        else:
            collection5.append("null")
    temp_df["Expressions Scores (2nd Lowest) - Tissue Type"] = collection5

    collection6 = []
    for ind in temp_df.index:
        var = temp_df["Expressions Scores (3rd Lowest)"][ind]
        if var.find("See") != -1:
            x = re.findall("\d+\.\d+", var)
            start = var.find(x[0])
            y = var[:start]
            collection6.append(y)
            temp_df["Expressions Scores (3rd Lowest)"][ind] = x[0]
        else:
            collection6.append("null")
    temp_df["Expressions Scores (3rd Lowest) - Tissue Type"] = collection6

    temp_df = temp_df[list(("Protein Name", "Signal Sequence", "Expressions Scores (Highest)",
                             "Expressions Scores (Highest) - Tissue Type", "Expressions Scores (2nd Highest)",
                             "Expressions Scores (2nd Highest) - Tissue Type", "Expressions Scores (3rd Highest)",
                             "Expressions Scores (3rd Highest) - Tissue Type", "Expressions Scores (Lowest)",
                             "Expressions Scores (Lowest) - Tissue Type", "Expressions Scores (2nd Lowest)",
                             "Expressions Scores (2nd Lowest) - Tissue Type", "Expressions Scores (3rd Lowest)",
                             "Expressions Scores (3rd Lowest) - Tissue Type", "UniPROT ID"))]

#Adding protein sequences and expression scores to respective dictionaries based on subcellular location
def compile(l1, l2, type):
    time.sleep(0.5)
    print(type)
    time.sleep(0.5)
    if type == "Membrane":
        print(name)
        mem_dict["Protein Name"].append(name)
        # print(p_id)
        # mem_dict["UniPROT ID"].append(p_id)
        print(sequence)
        mem_dict["Signal Sequence"].append(sequence)
        print(l1)
        print(l2)
        #Use "checker" variables to prevent "out of index" error
        checker1 = len(l1)
        checker2 = len(l2)
        print(checker1)
        print(checker2)
        if 1 <= checker1:
            mem_dict["Expressions Scores (Highest)"].append(l1[0])
        else:
            mem_dict["Expressions Scores (Highest)"].append("null")
        if 2 <= checker1:
            mem_dict["Expressions Scores (2nd Highest)"].append(l1[1])
        else:
            mem_dict["Expressions Scores (2nd Highest)"].append("null")
        if 3 <= checker1:
            mem_dict["Expressions Scores (3rd Highest)"].append(l1[2])
        else:
            mem_dict["Expressions Scores (3rd Highest)"].append("null")
        if 1 <= checker2:
            mem_dict["Expressions Scores (Lowest)"].append(l2[0])
        else:
            mem_dict["Expressions Scores (Lowest)"].append("null")
        if 2 <= checker2:
            mem_dict["Expressions Scores (2nd Lowest)"].append(l2[1])
        else:
            mem_dict["Expressions Scores (2nd Lowest)"].append("null")
        if 3 <= checker2:
            mem_dict["Expressions Scores (3rd Lowest)"].append(l2[2])
        else:
            mem_dict["Expressions Scores (3rd Lowest)"].append("null")
        print(mem_dict)
    if type == "Surface":
        print(name)
        surf_dict["Protein Name"].append(name)
        print(p_id2)
        surf_dict["UniPROT ID"].append(p_id2)
        print(sequence)
        surf_dict["Signal Sequence"].append(sequence)
        print(l1)
        print(l2)
        #Use "checker" variables to prevent "out of index" error
        checker1 = len(l1)
        checker2 = len(l2)
        print(checker1)
        print(checker2)
        if 1 <= checker1:
            surf_dict["Expressions Scores (Highest)"].append(l1[0])
        else:
            surf_dict["Expressions Scores (Highest)"].append("null")
        if 2 <= checker1:
            surf_dict["Expressions Scores (2nd Highest)"].append(l1[1])
        else:
            surf_dict["Expressions Scores (2nd Highest)"].append("null")
        if 3 <= checker1:
            surf_dict["Expressions Scores (3rd Highest)"].append(l1[2])
        else:
            surf_dict["Expressions Scores (3rd Highest)"].append("null")
        if 1 <= checker2:
            surf_dict["Expressions Scores (Lowest)"].append(l2[0])
        else:
            surf_dict["Expressions Scores (Lowest)"].append("null")
        if 2 <= checker2:
            surf_dict["Expressions Scores (2nd Lowest)"].append(l2[1])
        else:
            surf_dict["Expressions Scores (2nd Lowest)"].append("null")
        if 3 <= checker2:
            surf_dict["Expressions Scores (3rd Lowest)"].append(l2[2])
        else:
            surf_dict["Expressions Scores (3rd Lowest)"].append("null")
        print(surf_dict)
    if type == "ECM":
        print(name)
        ecm_dict["Protein Name"].append(name)
        print(p_id3)
        ecm_dict["UniPROT ID"].append(p_id3)
        print(sequence)
        ecm_dict["Signal Sequence"].append(sequence)
        print(l1)
        print(l2)
        #Use "checker" variables to prevent "out of index" error
        checker1 = len(l1)
        checker2 = len(l2)
        print(checker1)
        print(checker2)
        if 1 <= checker1:
            ecm_dict["Expressions Scores (Highest)"].append(l1[0])
        else:
            ecm_dict["Expressions Scores (Highest)"].append("null")
        if 2 <= checker1:
            ecm_dict["Expressions Scores (2nd Highest)"].append(l1[1])
        else:
            ecm_dict["Expressions Scores (2nd Highest)"].append("null")
        if 3 <= checker1:
            ecm_dict["Expressions Scores (3rd Highest)"].append(l1[2])
        else:
            ecm_dict["Expressions Scores (3rd Highest)"].append("null")
        if 1 <= checker2:
            ecm_dict["Expressions Scores (Lowest)"].append(l2[0])
        else:
            ecm_dict["Expressions Scores (Lowest)"].append("null")
        if 2 <= checker2:
            ecm_dict["Expressions Scores (2nd Lowest)"].append(l2[1])
        else:
            ecm_dict["Expressions Scores (2nd Lowest)"].append("null")
        if 3 <= checker2:
            ecm_dict["Expressions Scores (3rd Lowest)"].append(l2[2])
        else:
            ecm_dict["Expressions Scores (3rd Lowest)"].append("null")
        print(ecm_dict)

#Enter BGEE database to access expression scores for the given protein
def bgee_search(bgee_id, type):
    type = type
    print("Bgee start")
    bgee_url = "https://www.bgee.org/gene/" + bgee_id + "/"
    driver.get(bgee_url)
    time.sleep(5)
    print(driver.page_source.encode("utf-8"))
    #element1 contains html code of the highest expression scores and those respective tissue types
    element1 = driver.find_elements(By.XPATH, '//*[@id="gene-body"]/div[3]/div[5]/table')
    lines1 = []
    lines2 = []
    for value1 in element1:
        temp1 = value1.text
        count=0
        for line in temp1.splitlines():
            count+=1
            if (count > 4) and (count < 8):
                print(count)
                print(line)
                lines1.append(line)

    #element2 contains html code of the lowest expression scores and those respective tissue types
    element2 = driver.find_elements(By.XPATH, '//*[@id="gene-body"]/div[4]/div[5]/table')
    for value2 in element2:
        temp2 = value2.text
        count=0
        for line in temp2.splitlines():
            count += 1
            if (count > 4) and (count < 8):
                print(count)
                print(line)
                lines2.append(line)
    compile(lines1, lines2, type)

    print("Bgee done")
    time.sleep(0.5)

#Enter the UniProt profile for the given protein to access official name, signal peptide sequence, and BGEE ID for next steps
def uniprot_search(prot_id, type):
    type = type
    prot_url = "https://uniprot.org/uniprotkb/" + prot_id + "/entry"
    driver.get(prot_url)
    print("Entered protein site")
    time.sleep(2)
    driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")
    print("Scrolled, protein source code incoming")
    time.sleep(2)
    raw_html = driver.page_source.encode("utf-8")
    print(raw_html)

    time.sleep(0.5)
    print("Searching for official name")
    prot_name = driver.find_elements(By.XPATH, '//*[@id="root"]/div/div/div/main/ul/li[1]/div/div[2]/strong')
    for pn in prot_name:
        global name
        name = pn.text
        print(name)

    time.sleep(0.5)
    print("Searching for Signal Peptide Sequence")
    signal_peptide = driver.find_elements(By.XPATH, '//*[@id="ptm_processing"]/div/div[2]/protvista-manager/div[2]/protvista-datatable/table/tbody/tr[2]/td')
    for sp in signal_peptide:
        global sequence
        sequence = sp.get_attribute("innerText")
        print(sequence)

    time.sleep(0.5)
    print("Searching for BGEE link")
    bgee_list = driver.find_elements(By.XPATH, '//*[@id="expression"]/div/div[2]/ul[1]/li[1]/div/div[2]/ul/li/a')
    for bg in bgee_list:
        bg_id = bg.text
        print(bg_id)
        bgee_search(bg_id, type)

#Initialize Selenium webdriver for web scraping and data gathering purposes
chrome_options = Options()
chrome_options.add_argument("--headless")
chrome_options.add_argument('--no-sandbox')
chrome_options.add_argument('--disable-dev-shm-usage')
chrome_options.add_argument('enable-logging')
driver = webdriver.Chrome(options=chrome_options)

#Read CSV files of downloaded UniProt collections containing cell surface (surf_ex), ECM (ecm_ex), and cell membrane (mem_ex) proteins
surf_ex = pd.read_excel("C:/user/folder/surface_proteins.xlsx") #replace with appropriate filepath
ecm_ex = pd.read_excel("C:/user/folder/ECM_proteins.xlsx") #replace with appropriate filepath
mem_ex = pd.read_excel("C:/user/folder/membrane_proteins.xlsx") #replace with appropriate filepath

#Iterate through membrane proteins, using the unique protein ID to pass into uniprot_search function
i=0
for ind in mem_ex.index:
    i += 1
    p_name = mem_ex["Protein names"][ind]
    p_id = mem_ex["Entry"][ind]
    uniprot_search(p_id, "Membrane")
    print(i)

time.sleep(2)

#Iterate through surface proteins, using the unique protein ID to pass into uniprot_search function
j=0
global p_id2
for ind in surf_ex.index:
    j += 1
    p_name = surf_ex["Protein names"][ind]
    p_id2 = surf_ex["Entry"][ind]
    uniprot_search(p_id2, "Surface")
    print(j)

time.sleep(2)

#Iterate through ECM proteins, using the unique protein ID to pass into uniprot_search function
k=0
global p_id3
for ind in ecm_ex.index:
    k += 1
    p_name = ecm_ex["Protein names"][ind]
    p_id3 = ecm_ex["Entry"][ind]
    uniprot_search(p_id3, "ECM")
    print(k)
        
time.sleep(2)

driver.quit()
print("All Done")

#Create dataframes using each dictionary (mem_dict, surf_dict, ecm_dict)
df = pd.DataFrame(data=mem_dict)
print(df)
#Send database containing membrane proteins to be filtered and reorganized through filter_sort function
filter_sort(df)
df["Location"] = "Cell Membrane"
df = df[list(("Protein Name", "Signal Sequence", "Expressions Scores (Highest)",
                             "Expressions Scores (Highest) - Tissue Type", "Expressions Scores (2nd Highest)",
                             "Expressions Scores (2nd Highest) - Tissue Type", "Expressions Scores (3rd Highest)",
                             "Expressions Scores (3rd Highest) - Tissue Type", "Expressions Scores (Lowest)",
                             "Expressions Scores (Lowest) - Tissue Type", "Expressions Scores (2nd Lowest)",
                             "Expressions Scores (2nd Lowest) - Tissue Type", "Expressions Scores (3rd Lowest)",
                             "Expressions Scores (3rd Lowest) - Tissue Type", "Location", "UniPROT ID"))]

df2 = pd.DataFrame(data=surf_dict)
print(df2)
#Send database containing surface proteins to be filtered and reorganized through filter_sort function
filter_sort(df2)
df2["Location"] = "Cell Surface"
df2 = df2[list(("Protein Name", "Signal Sequence", "Expressions Scores (Highest)",
                             "Expressions Scores (Highest) - Tissue Type", "Expressions Scores (2nd Highest)",
                             "Expressions Scores (2nd Highest) - Tissue Type", "Expressions Scores (3rd Highest)",
                             "Expressions Scores (3rd Highest) - Tissue Type", "Expressions Scores (Lowest)",
                             "Expressions Scores (Lowest) - Tissue Type", "Expressions Scores (2nd Lowest)",
                             "Expressions Scores (2nd Lowest) - Tissue Type", "Expressions Scores (3rd Lowest)",
                             "Expressions Scores (3rd Lowest) - Tissue Type", "Location", "UniPROT ID"))]


df3 = pd.DataFrame(data=ecm_dict)
print(df3)
#Send database containing ECM proteins to be filtered and reorganized through filter_sort function
filter_sort(df3)
df3["Location"] = "Extracellular Matrix"
df3 = df3[list(("Protein Name", "Signal Sequence", "Expressions Scores (Highest)",
                             "Expressions Scores (Highest) - Tissue Type", "Expressions Scores (2nd Highest)",
                             "Expressions Scores (2nd Highest) - Tissue Type", "Expressions Scores (3rd Highest)",
                             "Expressions Scores (3rd Highest) - Tissue Type", "Expressions Scores (Lowest)",
                             "Expressions Scores (Lowest) - Tissue Type", "Expressions Scores (2nd Lowest)",
                             "Expressions Scores (2nd Lowest) - Tissue Type", "Expressions Scores (3rd Lowest)",
                             "Expressions Scores (3rd Lowest) - Tissue Type", "Location", "UniPROT ID"))]


#Print the three databases and the concatenated database to confirm successful output
print(db)
print(df2)
print(df3)

result = pd.concat([db, df2, df3], axis=0, ignore_index=True)
print(result)

result.to_csv("C:/user/folder/final_database.csv") #replace with appropriate filepath
