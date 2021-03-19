from urllib.request import urlopen
from urllib.error import HTTPError
import time
import pandas as pd
try:
    import pybel
except:
    from openbabel import pybel

# to use unverified ssl you can add this to your code
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

def populate_molidentifiers(df):
    """ Function to fetch the following molecule identifiers for input molecules from public sources (currently uses only https://cactus.nci.nih.gov/chemical/structure?identifier=&representation=names):

        1) SMILES (Simplified Molecular Input Line Entry System) notation
        2) Standard InChI (International Chemical Identifier)
        3) Standard InChI Key (27 character hashed International Chemical Identifier)
        4) Chemical formula
        5) IUPAC (International Union of Pure and Applied Chemistry) name
        6) other names (all other names for the given compound, including common names)
        7) CAS (Chemical Abstracts Services) registry number

    Parameters
    ----------
    df: pandas dataframe
        Molecule data that needs to be displayed on the molecule page. 

    Returns
    -------
    df: pandas dataframe
        molecule identifiers that are not included in col_list for the input molecules
    """
    identifier_dict = {'SMILES':'smiles','Standard_Inchi_Key':'stdinchikey','Standard_Inchi':'stdinchi','Chemical_Formula':'formula','IUPAC_Name':'iupac_name','Other_name':'names','CAS_Registry_Number':'cas'}
    from_pybel = {'SMILES':'can','Standard_Inchi_Key':'inchikey','Standard_Inchi':'inchi'}
    added = []
    url = ''
    m = ''
    not_none = ''
    for mol in identifier_dict.keys():
        if df[mol][0] != 'NONE':
            if mol in from_pybel.keys():
                #ans = df[col_list[0]].apply(lambda x:urlopen(url+x+'/'+identifier_dict[mol]).read().decode('utf8'))
                m = pybel.readstring(from_pybel[mol],df[mol][0])
            else:
                url='http://cactus.nci.nih.gov/chemical/structure/'+df[mol][0]+'/'
            not_none=mol
            break
    
    for mol in identifier_dict.keys():
        if df[mol][0] == 'NONE':
            if mol in from_pybel.keys():
                try:
                    added.append(mol)
                    ans = m.write(from_pybel[mol]).strip()
                    df[mol][0] = ans
                except TypeError:
                    continue
            else:
                if url == '':
                    new_url = 'http://cactus.nci.nih.gov/chemical/structure/'+df[not_none][0]+'/' + identifier_dict[mol]
                else:
                    new_url = url + identifier_dict[mol]
                try:
                    added.append(mol)
                    print(added)
                    print(url)
                    ans = urlopen(new_url).read().decode('utf8')
                    df[mol][0] = ans
                except HTTPError:
                    continue
                    #print('The {} for the {}({}) is either unavailable or \nthe value of {} provided by you is incorrect'.format(mol,col_list[0],row,row))
    #mol_tup = []
    #for i in range(len(df)):
    #    single_mol = []
    #    for k in molidentifiers.keys():
    #        single_mol.append(molidentifiers[k][i])
    #    mol_tup.append(single_mol)
    
    #print(mol_tup)
    return df,list(set(added).intersection(set(df.columns)))

                
#if 'inchi' in mol:
    #ans=ans[ans.index('=')+1:]