import sys
import filedownload
import requests
import os
import logging
import json
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.common.by import By

logger = logging.getLogger('main.InterEvol')

'''
This file contains the function that requests an alignment between two
PDB files in InterEvolAlign.
The central function is runAlign and is still in progress : it is not fully
automatic and for the moment, there is no way for the user to tinker the
parameters of the alignment or extract other informations that the PDB files
of the resulting interologs.
For the moment, the __main__ requires three arguments :
    - dir1 : the directory containing the PDB files of the subunits
    - dir2 : the directory containing the FASTA sequences
    - dir3 : the directory containing the PDB files of the interologs
    - receptor : the filename of the receptor
    - ligand : the filename of the ligand

!! CAUTION !!
For the moment, the runAlign function is not temporized to not send a request
to the InterEvol website will one is already running, and it could be possible to send
with some modifications many requests at once to the website : do not do that or you'll risk IP-ban
from the website. Only use runAlign do to a request at once.
!! CAUTION !!
'''
    
def runAlign(file1,file2,interDirectory):
    """
    Send a alignment request to the InterEvol database
    input : two FASTA files
    output : a dictionnary with the keys corresponding to the interologs
    and the values corresponding to the sequence identities of the ligand
    and the receptor
    """
    
    url = "http://biodev.cea.fr/interevol/interevalign.aspx"
    
    """
    Firstly, writes the content of the FASTA file in two variables fasta1 and fasta2
    """
    with open(file1, 'r') as content_file:
        fasta1 = content_file.read()
    with open(file2, 'r') as content_file:
        fasta2 = content_file.read()
    
    """
    Creates a webdriver object named browser to open the website
    Notes that with part of the script requires the installation of geckodriver
    Latest releases of geckodriver : https://github.com/mozilla/geckodriver/releases
    On Unix system, extract the executable to /usr/bin
    Then make geckodriver executable with chmod +x geckodriver
    If you run into an error with geckodriver, an easy fix is to retrograde your selenium version
    For that, simply use the command "pip install selenium==2.53.6
    """

    browser = webdriver.Firefox()
    browser.implicitly_wait(10)
    browser.get(url)
    
    """
    Fills the HTML forms of interests, mainly :
        - 'seqReceptor' textarea with the FASTA sequence contained in fasta1
        - 'seqLigand' textarea with the FASTA sequence contained in fasta2
        - Submit the job by clicking on the "Submit" button
    """
    browser.find_element_by_id('seqReceptor').send_keys(fasta1) 
    browser.find_element_by_id('seqLigand').send_keys(fasta2)
    browser.find_element_by_id('btnRun').click() 
    
    
    delay = 3600 #in seconds
    logger.debug("Starting alignment")
    
    try:
        WebDriverWait(browser, delay).until(EC.title_is('InterEvolAlign results page'))
        #The browser will wait the delay or until the results page appears before anymore script is run
        logger.debug("Page is ready!") #If the page appears, the rest of the script is run
    except TimeoutException:
        logger.error("Loading took too much time!")
        browser.quit()
        return #Otherwise, the function is stopped
        
    """
    Looks for a href tag in the source code with "interid" (interology)
    """
    interolog = "interid"  
    InterologList = browser.find_elements_by_xpath('//a[contains(@href, "%s")]' % interolog)        
    PDBid=list()
    rchain=list()
    lchain=list()
    result = {}
    
    """
    Extract the four letter PDB id of interelog domains
    """
    for element in InterologList:
        e=element.text
        r = e[5:-1]
        l = e[6:]
        logger.debug('{} {} {}'.format(e,r,l))
        PDBid.append(e)
        rchain.append(r)
        lchain.append(l)
    
    """
    First step is to find the table containing the interologs
    """
    table_id = browser.find_element_by_id('GridView3')
    rows = table_id.find_elements_by_tag_name("tr") # get all of the rows in the table
    rows.pop(0)
    
    """
    The next step is to extract the identity percentages    
    """
    Idllist=list()
    Idrlist=list()
    
    for row in rows:
        #Get the columns      
        coll = row.find_elements_by_tag_name("td")[2] #note: index start from 0, 1 is col 2
        colr = row.find_elements_by_tag_name("td")[3]
        #logger.debugs text from the element
        logger.debug((coll.text))
        logger.debug((colr.text))
        Idllist.append(coll.text)
        Idrlist.append(colr.text)
    
    """
    The interologs and the identity percentages are stored in a dictionnary
    The keys correspond to the interologs PDB ID
    The values of each keys is a list containing the identity percentages
    If the identity with of the receptor and the ligand is 100%, it means that
    the complex already exists
    Therefore, the interolog is not downloaded and the PDB id is removed from the PDBid list
    The values also contains the homologous chain of the receptor
    and the ligand on the interolog complex
    """
    
    j=-1
    for i in range(0,len(Idllist)):
        j=j+1
        if((Idllist[i] == "100%") and (Idrlist[i] == "100%")):
            logger.info("The complex seems to have already been caracterized. (ref {})".format(PDBid[i]))
            del PDBid[j]
            j=j-1
        else:
            result[PDBid[j]]={'idn_l':Idllist[i],'idn_r':Idrlist[i],'chn_l':lchain[i],'chn_r':rchain[i]}
        
    logger.debug(PDBid)
    logger.debug(result)

    """
    Downloads the PDB file of interologs from InterEvol
    """
    logger.debug("Downloading PDF file of each interolog.")
    for element in PDBid:
        url = "http://biodev.cea.fr/interevol/interevol.aspx?interid="+element+"#footer"
        browser.get(url)
        WebDriverWait(browser, 10).until(EC.text_to_be_present_in_element((By.ID, "downloadAlignment"), element[:4]))
        logger.debug(browser.current_url)
        tag = browser.find_element_by_id("downloadAlignment")
        WebDriverWait(browser, 10).until(EC.presence_of_element_located((By.ID, "downloadAlignment")))
        url2 = tag.get_attribute("href")
        logger.debug(url2)
        r = requests.get(url2)
        with open(os.path.join(interDirectory,element+".pdb"), "wb") as code:
            code.write(r.content)
        
        
    logger.debug("Finished downloading")
    logger.debug((browser.current_url))
    logger.info("{:02d} interologs are found. ({})".format(len(PDBid),' '.join(PDBid)))
    return(result)
        
if __name__ == '__main__':
    print("Begin")
    try:
        #The directory containing the PDB files of the subunits
        subunitsDirectory = sys.argv[sys.argv.index("-dir1")+1]
        #The directory where the FASTA sequences will be stored
        fastaDirectory = sys.argv[sys.argv.index("-dir2")+1]
        #The directory where the PDB files of the interologs will be stored
        interDirectory = sys.argv[sys.argv.index("-dir3")+1]
                
    except:    
        print("ERROR: please enter the names of the directories\n")
        sys.exit()
        
    try:
        #The name of the PDB file of the receptor
        ligand = sys.argv[sys.argv.index("-receptor")+1]
        #The name of the PDB file of the ligand
        receptor = sys.argv[sys.argv.index("-ligand")+1]
    except:
        print("ERROR : specified file does not exist\n")
        sys.exit()
        
    
    print("Creating FASTA sequence")
    ligf = filedownload.FASTAfromPDB(ligand, subunitsDirectory, fastaDirectory)
    recf = filedownload.FASTAfromPDB(receptor, subunitsDirectory, fastaDirectory)
    
    print("Alignment function")
    
    """
    Keys of dictionnary : ensemble of interologs found by InterEvol (format example : 2za4_AB)
    Values of dictionnary : a list containing in order 
		1) The sequence identity with the receptor in percentage
		2) The sequence identity with the ligand in percentage
		3) The homologous chain on the interolog corresponding to the receptor given in input
	    4) The homologous chain on the interolog corresponding to the ligand given in input
	"""
    dico = runAlign(recf,ligf,interDirectory)
    
