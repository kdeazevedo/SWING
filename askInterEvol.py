import sys
import filedownload
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException


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

!! CAUTION !!
For the moment, the runAlign function is not temporized to not send a request
to the InterEvol website will one is already running, and it is possible to send
many requests at once to the website : do not do that or you'll risk IP-ban
from the website. Only use runAlign do to a request at once.
!! CAUTION !!
'''
    
def runAlign(file1,file2):
    """
    Send a alignment request to the InterEvol database
    input : two FASTA files
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
    browser.implicitly_wait(5)
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
    
    
    delay = 600 #in seconds
    print "Starting alignment"
    
    try:
        WebDriverWait(browser, delay).until(EC.title_is('InterEvolAlign results page'))
        #The browser will wait the delay or until the results page appears before anymore script is run
        print "Page is ready!" #If the page appears, the rest of the script is run
    except TimeoutException:
        print "Loading took too much time!"
        browser.quit()
        return #Otherwise, the function is stopped
    
    browser.implicitly_wait(5)
    
    """
    Looks for a href tag in the source code with "interid" (interology)
    """
    interolog = "interid"  
    InterologList = browser.find_elements_by_xpath('//a[contains(@href, "%s")]' % interolog)        
    PDBid=list()
    
    """
    Extract the four letter PDB id of interelog domains
    """
    for element in InterologList:
        e=element.text
        e=e[:4]
        PDBid.append(e)
        print(e)

    print(browser.current_url)
    return(PDBid)
        
if __name__ == '__main__':
    try:
        #The directory containing the PDB files of the subunits
        subunitsDirectory = sys.argv[sys.argv.index("-dir1")+1]
        #The directory containing the FASTA sequences of the subunits
        fastaDirectory = sys.argv[sys.argv.index("-dir2")+1]
        #The directory where the PDB files of the interologs will be stored
        interDirectory = sys.argv[sys.argv.index("-dir3")+1]
        
    except:    
        print "ERROR: please enter the names of the directories\n"
        sys.exit()
    
    filedownload.parseDirectory(subunitsDirectory, fastaDirectory)
    
    for files in fastaDirectory:
        '''
        For the moment, we'll use only these two fasta sequences for testing purposes
        The script will be modified to not have to manually enter the files names
        '''
        file1=fastaDirectory+"2KF4.fasta.txt"
        file2=fastaDirectory+"1BTA.fasta.txt"
    
    PDBid = runAlign(file1,file2)
    filedownload.downloadPDB(PDBid, interDirectory)
    

