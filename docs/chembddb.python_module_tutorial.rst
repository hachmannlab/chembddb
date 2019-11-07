ChemBDDB Python Module Tutorial
===============================

The setup, insert, and delete modules of chembddb can also be called using python code, making it easier to use and simple. 

Step 1: import the chembddb module and initialize variables for mysql credentials
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.. code:: python

    import chembddb
    host = '127.0.0.1'
    username = 'root'
    password = 'yourpassword'

Step 2: call the setup module
++++++++++++++++++++++++++++++
.. code:: python

    # uncomment the following line if you want to know more about the setup module
    # print(help(chembddb.setup))

    chembddb.setup(host=host, user=username, pw=password, db='abc')

Step 3.1: call the insert module and use the path to the configuration and data files
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.. code:: python

    # uncomment the following line if you want to know more about the setup module
    # print(help(chembddb.insert))

    chembddb.insert(host=host, user=username, pw=password, db='abc',smi_col='smiles',mol_identifier='name',conf_file='../test-files/config_hsp.csv',data_file='../test-files/benzene.csv')

Step 3.2: call the insert module and use pandas dataframes for the configuration and data files
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.. code:: python

    import pandas as pd
    conf = pd.read_csv('../test-files/config_hsp.csv')
    data = pd.read_csv('../test-files/benzene.csv')

    # uncomment the following line if you want to know more about the setup module
    # print(help(chembddb.insert))

    chembddb.insert(host=host, user=username, pw=password, db='abc',smi_col='smiles',mol_identifier='name',conf_file=conf,data_file=data)

Step 4: call the delete module
+++++++++++++++++++++++++++++++
.. code:: python

    # uncomment the following line if you want to know more about the setup module
    # print(help(chembddb.delete))
    
    chembddb.delete(host=host,user=username,pw=password,db='abc')
