from flask import Flask, render_template, url_for, request,redirect
import pymysql
import os
import sys
import pandas as pd
from copy import deepcopy
from chembddb import molidentfiers
from chembddb.units import insert_unit_list,fetch_unit_list, create_unit_list, unit_converter
import json

try:
    import pybel
except:
    from openbabel import pybel

from flask import send_from_directory
import numpy as np
import time

#cur = ''
#con = ''
all_dbs = []
app = Flask(__name__)
upload_directory=os.getcwd()
app.config['UPLOAD FOLDER']=upload_directory

def post_process(sql, from_db):
    """
    Helper function to post process the results from the database

    Parameters
    ----------
    sql: str
        sql query that was created by the search function
    from_db: tuple of tuples
        results from the database

    Returns
    -------
    data: pandas dataframe
        dataframe containing the processed results ready to be displayed
    columns: list of str
        list of formatted and cleaned column names
    """
    abc = 0
    if 'MW' in sql and 'property' not in sql:
        data = pd.DataFrame(list(from_db), columns=['ID','SMILES','MW'])
        columns = list(data.columns)
    elif 'MW' not in sql and 'property' not in sql:
        data = pd.DataFrame(list(from_db),columns = ['ID','SMILES'])
        columns = list(data.columns)
    else:
        data = pd.DataFrame(list(from_db), columns=['Molecule_id','SMILES','Method','Functional','Basis_set','forcefield','Property','Value'])
        data['ID_SMI']=data['Molecule_id'].astype(str)+','+data['SMILES']
        data['Property']=data['Property']+'-' +data['Method']+'('+data['Functional']+'/'+data['Basis_set']+')('+data['forcefield']+')'
        data = data[data.columns[-3:]]
        data=data.pivot_table(index='ID_SMI',columns='Property',values='Value')
        data = data.reset_index()
        data[['ID','SMILES']]=data['ID_SMI'].str.split(',',expand=True)
        columns=['ID','SMILES']
        for i in data.columns[1:-2]:
            columns.append(i)
        data=data[columns]
        columns=[c.replace('(NA/NA)','') for c in columns]
        columns=[c.replace('(na/na)','') for c in columns]
        columns=[c.replace('(NA)','') for c in columns]
        columns=[c.replace('(na)','') for c in columns]
    return data, columns

@app.route('/')
def begin():
    return redirect(url_for('connect'))

def connect_mysql(host,user,pw):
    global cur, con
    try:
        con = pymysql.connect(host = host, user=user, password = pw)
        cur = con.cursor()
        cur.execute('show databases;')
        all_dbs_tup=cur.fetchall()
        all_dbs = []
        for i in all_dbs_tup:
            if '_chembddb' in i[0]:
                m=i[0]
                all_dbs.append((m[:-9],))
        return cur,all_dbs
    except:
        return 'invalid','credentials'
    
@app.route('/connect',methods=['GET','POST'])
def connect():
    """establishes mysql connection based on the credentials provided during setup

    Parameters
    ----------
    Returns
    -------
    cur: cursor object
        pointer to the MySQL database; used to execute all SQL queries
    """ 
    global cur,all_dbs, unit_list,con
    if request.method=='POST':
        cred=request.form
        cred = cred.to_dict(flat=False)
        cur,all_dbs = connect_mysql(host = cred['host'][0], user=cred['username'][0], pw = cred['password'][0])
        if cur == 'invalid' and all_dbs == 'credentials':
            return render_template('connect.html',err_msg='Invalid Credentials. Did not connect to MySQL')
        else:
            print(all_dbs)
            if any('unit_list' in i for i in all_dbs):
                all_dbs.pop(all_dbs.index(('unit_list',)))
                unit_list = fetch_unit_list(cur)
            else:
                unit_list = create_unit_list(cur,con)
            return render_template('connect.html',success_msg='Connection Established',host=cred['host'][0],user=cred['username'][0],password=cred['password'][0],all_dbs=all_dbs)
    else:
        return render_template('connect.html')
    
@app.route('/setup',methods=['GET','POST'])
def setup(host=-1,user='',pw='',db=''):
    """
    Function to setup the database with the chembddb schema

    Parameters
    ----------
    host: str default=''
        the hostname is the domain name or server name
    user: str default=''
        the username for mysql
    pw: str default=''
        the password for mysql
    db: str default=''
        the name of the database that needs to be set up

    """
    if host != -1:
        # for python module
        b, a = connect_mysql(host=host,user=user,pw=pw)
        if b == 'invalid' and a == 'credentials':
            return 'invalid credentials'
        else:
            db = db +'_chembddb'
    elif request.method=='POST':
        # for UI
        db_details=request.form
        db_details=db_details.to_dict(flat=False)
        db=db_details['dbname'][0]+'_chembddb'
    else:
        # Default landing page for setup
        all_dbs=[]
        cur.execute('show databases;')
        all_dbs_tup=cur.fetchall()
        for i in all_dbs_tup:
            if '_chembddb' in i[0] and 'unit_list' not in i[0]:
                m=i[0]
                all_dbs.append((m[:-9],))

        return render_template('setup.html',all_dbs=all_dbs)

    all_dbs=[]
    cur.execute('show databases;')
    all_dbs_tup=cur.fetchall()
    for i in all_dbs_tup:
        if '_chembddb' in i[0] and 'unit_list' not in i[0]:
            m=i[0]
            all_dbs.append((m[:-9],))
    cur.execute('USE INFORMATION_SCHEMA')
    result=cur.execute('SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA WHERE SCHEMA_NAME=\'%s\''%db)
    if result == 0:
        cur.execute('CREATE DATABASE %s;'%db)
        cur.execute('USE %s;'%db)
        cur.execute('CREATE TABLE `%s`.`Molecule` (`id` INT NOT NULL AUTO_INCREMENT,`SMILES` VARCHAR(300) DEFAULT \'NONE\', `Standard_InChI` VARCHAR(400) DEFAULT \'NONE\',`Standard_InChI_Key` VARCHAR(100) DEFAULT \'NONE\',`CAS_Registry_Number` VARCHAR(200) DEFAULT \'NONE\',`IUPAC_Name` VARCHAR(400) DEFAULT \'NONE\',`Other_name` VARCHAR(1000) DEFAULT \'NONE\',`Chemical_Formula` VARCHAR(100) DEFAULT \'NONE\',`MW` FLOAT, PRIMARY KEY (`id`));'%db)
        # cur.execute('CREATE TABLE `%s`.`Credit`(`id` INT NOT NULL AUTO_INCREMENT,`DOI` VARCHAR(100) UNIQUE DEFAULT \'None\',`details` VARCHAR(100) DEFAULT \'None\',PRIMARY KEY (`id`));'%db)
        cur.execute('CREATE TABLE `%s`.`Property`(`id` INT NOT NULL AUTO_INCREMENT,`Property_str` VARCHAR(100) NOT NULL UNIQUE,`Unit` VARCHAR(100) NOT NULL,PRIMARY KEY (`id`));'%db)
        cur.execute('CREATE TABLE `%s`.`Model`(`id` INT NOT NULL AUTO_INCREMENT,`method_name` VARCHAR(100) NOT NULL UNIQUE,`options` VARCHAR(500),PRIMARY KEY (`id`));'%db)
        cur.execute('CREATE TABLE `%s`.`Functional`(`id` INT NOT NULL AUTO_INCREMENT,`name` VARCHAR(100) DEFAULT \'NONE\',PRIMARY KEY (`id`));'%db)
        cur.execute('CREATE TABLE `%s`.`Basis_set`(`id` INT NOT NULL AUTO_INCREMENT,`name` VARCHAR(100) DEFAULT \'NONE\',PRIMARY KEY (`id`));'%db)
        cur.execute('CREATE TABLE `%s`.`Forcefield`(`id` INT NOT NULL AUTO_INCREMENT,`name` VARCHAR(100) DEFAULT \'NONE\',PRIMARY KEY (`id`));'%db)
        # cur.execute('CREATE TABLE `%s`.`Topology`(`id` INT NOT NULL AUTO_INCREMENT,`geometry` VARCHAR(100) NOT NULL,`symbols` VARCHAR(100),`method` VARCHAR(100),`steps` INT,PRIMARY KEY (`id`));'%db)
        cur.execute('CREATE TABLE `%s`.`Value`(`id` INT NOT NULL AUTO_INCREMENT,`num_value` FLOAT NOT NULL,`model_id` INT NOT NULL,`property_id` INT NOT NULL,`molecule_id` INT NOT NULL,`functional_id` INT, `basis_id` INT,`forcefield_id` INT,PRIMARY KEY (`id`));'%db)
        cur.execute('CREATE TABLE `%s`.`Configuration`(`id` INT DEFAULT 0,`conf` VARCHAR(200) DEFAULT \'NONE\',`unit_dict` VARCHAR(500) DEFAULT \'NONE\',PRIMARY KEY (`id`));'%db)
        cur.execute('ALTER TABLE `%s`.`Value` ADD CONSTRAINT `Value_fk0` FOREIGN KEY (`model_id`) REFERENCES `Model`(`id`) on DELETE CASCADE;'%db)
        cur.execute('ALTER TABLE `%s`.`Value` ADD CONSTRAINT `Value_fk1` FOREIGN KEY (`property_id`) REFERENCES `Property`(`id`) on DELETE CASCADE;'%db)
        cur.execute('ALTER TABLE `%s`.`Value` ADD CONSTRAINT `Value_fk2` FOREIGN KEY (`molecule_id`) REFERENCES `Molecule`(`id`) on DELETE CASCADE;'%db)
        cur.execute('ALTER TABLE `%s`.`Value` ADD CONSTRAINT `Value_fk3` FOREIGN KEY (`functional_id`) REFERENCES `Functional`(`id`) on DELETE CASCADE;'%db)
        cur.execute('ALTER TABLE `%s`.`Value` ADD CONSTRAINT `Value_fk4` FOREIGN KEY (`basis_id`) REFERENCES `Basis_set`(`id`) on DELETE CASCADE;'%db)
        cur.execute('ALTER TABLE `%s`.`Value` ADD CONSTRAINT `Value_fk5` FOREIGN KEY (`forcefield_id`) REFERENCES `Forcefield`(`id`) on DELETE CASCADE;'%db)
        cur.execute('show databases;')
        all_dbs_tup=cur.fetchall()
        all_dbs=[]
        for i in all_dbs_tup:
            if '_chembddb' in i[0] and 'unit_list' not in i[0]:
                m=i[0]
                all_dbs.append((m[:-9],))
        if host == -1:
            # successful creation for UI
            return render_template('setup.html',dbname=db,all_dbs=all_dbs,success_msg='The database has been created.')
        else:
            # successful creation for python module
            return 'Success'
    else:
        if host == -1:
            # error handling for UI
            return render_template('setup.html',dbname=db,all_dbs=all_dbs,err_msg='Database already exists.')
        else:
            # error handling for python module
            return 'Failed! Database already exists.'

@app.route('/temp_insert',methods=['GET','POST'])
def temp_insert():
    global all_dbs, cur, data,db, mol_ids, con, snapshot, props, prop_names, unit_list, prop_type
    mi_cols=[]
    cur.execute('show databases;') 
    all_dbs_tup=cur.fetchall()
    all_dbs=[]
    for i in all_dbs_tup:
        if '_chembddb' in i[0] and 'unit_list' not in i[0]:
            m=i[0]
            all_dbs.append((m[:-9],))
    #print(request.form)
    if request.method=='POST' and 'upload_data' in request.form:
        # for UI with only the data file
        config_options=request.form
        config_options=config_options.to_dict(flat=False)
        db=config_options['dbname'][0]
        db = db+'_chembddb'
        files = request.files
        #files = files.to_dict(flat=False)
        cur.execute('USE {};'.format(db))
        cur.execute('SELECT ID,conf from Configuration')
        confs = cur.fetchall()
        conf = False
        if len(confs) > 0:
            #confs = confs[0][1].split(',')
            #confs = [confs[:7],confs[7:]]
            conf = True
        present = []
        identifiers = ['SMILES','Standard_Inchi_Key','Standard_Inchi','Chemical_Formula','IUPAC_Name','Other_name','CAS_Registry_Number']
        # to check which identifiers have entries in the database against them. The default value is 'NONE' so COUNT(DISTINCT) will give us number of unique entries for a particular column. If they are greater than 1 then the column has atleast one entry
        snapshot = []
        # IF SOMEONE MAKES A MISTAKE AND ADDS ONLY PROPERTY AND META DATA AND FORGETS TO ADD THE MOLECULES, THIS WILL NOT WORK.
        for id in identifiers:
            cur.execute('SELECT COUNT(DISTINCT '+ id +') FROM MOLECULE;')
            counts = cur.fetchall()
            #print(counts)
            if counts[0][0] > 1:
                present.append(id)
            # TODO: if value is not default and count is 1 then append to present
        if present !=[]:
            snapshot = [', '.join(p for p in present)]
            snapshot = [snapshot[0].split(',')]
            # check which properties are present in the database

            cur.execute('SELECT Property_str, Unit from Property;')
            properties = cur.fetchall()
            properties = ', '.join(i[0]+'('+i[1]+')' for i in properties)
            properties.replace(', na','')
            properties.replace('na,','')
            properties.replace('na','')
            
            properties = properties.split(',')
            snapshot.append(properties)

            cur.execute('SELECT method_name from Model;')
            methods = cur.fetchall()
            methods = ', '.join(i[0] for i in methods)
            methods = methods.replace(', na','')
            methods = methods.replace('na,','')
            methods = methods.replace('na','')

            methods = methods.split(',')
            snapshot.append(methods)

            cur.execute('SELECT name FROM Functional;')
            functionals = cur.fetchall()
            functionals = ', '.join(i[0] for i in functionals)
            functionals = functionals.replace(', na','')
            functionals = functionals.replace('na,','')
            functionals = functionals.replace('na','')
            
            functionals = functionals.split(',')
            snapshot.append(functionals)

            cur.execute('SELECT name FROM Basis_set;')
            basis = cur.fetchall()
            basis = ', '.join(i[0] for i in basis)
            basis = basis.replace(', na','')
            basis = basis.replace('na,','')
            basis = basis.replace('na','')
            
            basis = basis.split(',')
            snapshot.append(basis)

            cur.execute('SELECT name FROM Forcefield')
            forcefield = cur.fetchall()
            forcefield = ', '.join(i[0] for i in forcefield)
            forcefield = forcefield.replace(', na','')
            forcefield = forcefield.replace('na,','')
            forcefield = forcefield.replace('na','')
            forcefield = forcefield.replace(' ','')

            forcefield = forcefield.split(',')
            snapshot.append(forcefield)
        else:
            snapshot = False
        #print(snapshot)
        data_file = files['data_file']
        print(data_file.filename)
        if data_file.filename.rsplit('.',1)[1]!='csv':
            db.replace('_chembddb','')
            db=db.replace('_',' ')
            return render_template('temp_insert.html',all_dbs = all_dbs,title=db,err_msg='No data file provided or incorrect file format. (csv required)')
        else:
            data = pd.read_csv(data_file)
            cols = []
            for i in data.columns:
                cols.append(i.replace(' ','_'))
            print(data.columns)
            print(i)
            data.columns = cols
            return render_template('temp_insert.html',all_dbs=all_dbs,data_validated=True,cols = list(data.columns),conf=conf,snapshot=snapshot,snapshot_cols = ['Molecule Identifiers','Properties (Units)','Methods','Functionals','Basis Sets','Forcefields'])

    elif request.method == 'POST' and ('config' in request.form or 'use-config' in request.form):
        #print(data.head())
        config_options = request.form
        config_options=config_options.to_dict(flat=False)
        print(config_options)
            
        # loop throught he CSV file, check if the smiles value is in the table, if yes, fetch the corresponding id, same goes for property, same goes for method, if it does not exist, fetch the last id and create a new entry
        # populating and property table
        molecule_identifiers_cols = [] 
        molecule_identifiers = []
        cols = list(data.columns)

        for key in config_options.keys():
            if 'hidden' in key:
                molecule_identifiers_cols.append(config_options[key][0])
                molecule_identifiers.append(key[len('hidden_'):-3])
            else:
                continue
        #print(cols)
        print(molecule_identifiers)
        remaining_cols = list(set(cols)-set(molecule_identifiers_cols))
        mol_ids = {}
        for i in range(len(molecule_identifiers)):
            mol_ids[molecule_identifiers[i]] = molecule_identifiers_cols[i]

        if 'use-config' in config_options:
            #df = pd.read_csv(app.config['UPLOAD FOLDER']+'/'+db[:-9]+'_config.csv')
            #print(df)
            cur.execute('USE {};'.format(db))
            cur.execute('SELECT * from Configuration')
            confs = cur.fetchall()
            confs = confs[0][1].split(',')
            confs = [confs[:6],confs[6:]]
            print(confs)
            properties = []
            units = []
            methods = []
            functional = []
            basis = []
            forcefield = []
            for i in confs:
                properties.append(i[0])
                units.append(i[1])
                methods.append(i[2])
                functional.append(i[3])
                basis.append(i[4])
                forcefield.append(i[5])
            return render_template('temp_insert.html',all_dbs=all_dbs,conf_flag=True, l=len(properties), snapshot=snapshot,snapshot_cols=['Molecule Identifiers','Properties (Units)','Methods','Functionals','Basis Sets','Forcefields'],properties=properties,props = remaining_cols,prop_length=len(remaining_cols),title=db, unit=units,methods = methods, functional = functional, basis = basis, forcefield = forcefield)
        else:
            return render_template('temp_insert.html',props =remaining_cols, prop_length = len(remaining_cols), all_dbs=all_dbs, title=db,snapshot=snapshot,snapshot_cols = ['Molecule Identifiers','Properties (Units)','Methods','Functionals','Basis Sets','Forcefields'])

    elif request.method == 'POST' and ('meta-data' in request.form or 'download-submit' in request.form):
        meta_data = request.form
        meta_data = meta_data.to_dict(flat=False)
        props = []
        print(meta_data)

        for k in meta_data.keys():
            if 'prop_id' in k and 'type' not in k:
                props.append(meta_data[k][0])

        print(props)
        try:
            prop_type = meta_data['hidden_prop_id_type']
        except:
            prop_type = meta_data['prop_type_0']
        prop_names = meta_data['2_prop']
        print(prop_names)
        l = unit_list
        units = []
        for n in range(len(prop_names)):
            if prop_names[n].lower() in l.keys():
                units.append(list(l[prop_names[n].lower()].keys()))
            elif prop_type[n] == 'sol':
                units.append(list(l['solubility parameters'].keys()))
            elif prop_type[n] == 'energy':
                units.append(list(l['energy'].keys()))
            elif prop_type[n] == 'ratio':
                units.append(list(l['ratio'].keys()))
            else:
                units.append([''])
        to_drop=[]
        #print(units)
        return render_template('temp_insert.html',properties=prop_names,prop_length = len(prop_names),units = True, unit_list=units,title=db,all_dbs=all_dbs,snapshot=snapshot,snapshot_cols=['Molecule Identifiers','Properties (Units)','Methods','Functionals','Basis Sets','Forcefields'])
    
    elif request.method == 'POST' and 'meta-unit' in request.form:
        final_md = request.form
        final_md = final_md.to_dict(flat=False)
        print(unit_list)
        print(final_md)
        print(db)
        print(snapshot)
        print(mol_ids)
        print(data.columns)
        print(prop_names)
        print(props)
        unit_flag=False
        to_drop=[]
        units_from_ui = final_md['unit_id_0']
        cur.execute('Use {};'.format(db))
        #for md in final_md.keys():
        #    if 'unit_id' in md:
        #        if '(default)' in final_md[md]:
        #            unit_list.append(final_md[md][:-len('(default)')])
        #        else:
        #            unit_list.append(fi)
        
        #for p in range(len(prop_names)):
        #    if prop_names[p]+'('+final_md['unit']
        ## fetching properties from the table and checking them against the entered properties from the insert page. If they are not already existing in the table, insert them. 
        #entered_list=[]
        #cur.execute("SELECT Property_str from Property")
        #properties = cur.fetchall()
        #print(properties)
        #for prop, units in zip(props,meta_data['2_unit']):
        #    if any(prop in i for i in properties) or prop in entered_list:
        #        pass
        #    else:
        #        entered_list.append(prop)
        #        # TODO: try insert if does not exist using the SQL command
        #        print(entered_list)
        #        cur.execute("INSERT INTO Property(Property_str,Unit) VALUES(%s,%s)",[prop,units])
        #print('entered property')
        # populating the model table
        #entered_list=[]
        #for method in meta_data['2_method']:
        #    if any(method in i for i in models) or method in entered_list:
        #        pass
        #    else:
        #        entered_list.append(method)
        #        cur.execute("INSERT INTO Model(Method_name) VALUES(%s)",[method])
        #entered_list = []
        #for i in range(len(prop_names)):
        #    if '(default)' in final_md['2_unit'][i]:
        #        print(final_md['2_unit'][:-10])
        #        prop_unit = prop_names[i] + '('final_md['2_unit'][:-10]')'
        #        print(prop_unit)
        #        if prop_unit not in snapshot[1] and prop_unit not in entered_list:
        #            cur.execute('INSERT INTO Property(Property_str, Unit) VALUES(%s,%s)',[prop_names[i],final_md['2_unit'][:-10]])
        #            entered_list.append(prop_unit)
        #    else:
        #        if prop_names[i] in unit_list.keys():
        #            prop_unit = prop_names[i]+'('unit_list[prop_names[i]][0]+')'
        #            if prop_unit not in snapshot[1] and prop_unit not in entered_list:
        #                cur.execute('INSERT INTO Property(Property_str, Unit) VALUE(%s,%s)',[prop_names[i],unit_list[prop_names[i]][0]])
        #                entered_list.append(prop)
        #        else:
        #            # new unit
        # IF ',' in the unit then the unit is new (not in the unit)
        if snapshot!=False:
            # DB is not empty (there are some properties in the db)
            entered_list = []
            p = [s.split('(')[0].strip() for s in snapshot[1]] # properties in the database
            u = [s[s.index('('):s.index(')')] for s in snapshot[1]] # units in the database corresponding to that property
            for i in range(len(prop_names)):
                if (prop_names[i],units_from_ui[i]) not in entered_list:
                    print(prop_names[i])
                    print(unit_list.keys())
                    if prop_names[i] not in unit_list.keys() and prop_type[i] not in ['sol','energy']: # If property is not in our list
                        # creating unit entry for that property to add to our list and make that unit the default
                        d = {}
                        d[units_from_ui[i]+' (default)'] = 1.0
                        unit_list[prop_names[i]] = d
                        unit_flag=True
                        cur.execute('INSERT INTO Property(Property_str, unit) VALUE(%s,%s);',[prop_names[i],units_from_ui[i]])
                        # TODO: add to unit db and convert values
                    else: # if we already have the property in our list (of units)
                        if prop_names[i] not in p: # if the property is not in the database
                            if '(default)' in units_from_ui[i]:
                                u1 = units_from_ui[i]
                                cur.execute('INSERT INTO Property(Property_str,Unit) VALUES(%s,%s)',[prop_names[i],u1[:u1.index('(')-1]])
                            else:
                                if prop_type[i] == 'sol':
                                    all_units = list(unit_list['solubility parameters'].keys())
                                    def_unit = all_units[0]
                                elif prop_type[i] == 'energy':
                                    all_units = list(unit_list['energy'].keys())
                                    def_units = all_units[0]
                                else:
                                    all_units = list(unit_list[prop_names[i]].keys())
                                    def_unit = all_units[0] # fetching the default unit

                                if units_from_ui[i] not in all_units: # if we do not have this unit listed for the given property in our list
                                    conv_factor = float(units_from_ui[i].split(',')[1])
                                    units_from_ui[i] = units_from_ui[i].split(',')[0]
                                    if prop_type[i] == 'sol':
                                        unit_list['solubility parameters'][units_from_ui[i]] = conv_factor
                                    elif prop_type[i] == 'energy':
                                        unit_list['energy'][units_from_ui[i]] = conv_factor
                                    else:
                                        unit_list[prop_names[i]][units_from_ui[i]] = conv_factor
                                    unit_flag=True
                                else:
                                    # TODO: convert
                                    if prop_type[i] == 'sol':
                                        conv_factor = unit_list['solubility parameters'][units_from_ui[i]]
                                    elif prop_type[i] == 'energy':
                                        conv_factor = unit_list['energy'][units_from_ui]
                                    else:
                                        conv_factor = unit_list[prop_names[i]][units_from_ui[i]]
                                        
                                data[props[i]] = data[props[i]] * conv_factor
                                cur.execute('INSERT INTO Property(Property_str,Unit) VALUES(%s,%s)',[prop_names[i],def_unit])
                                print(def_unit)
                        else: # If the property is already in the database
                            if '(default)' not in units_from_ui[i]:
                                def_unit = list(unit_list[prop_names[i]].keys())[0] 
                                if units_from_ui[i] not in unit_list[prop_names[i]].keys():
                                    unit_list[prop_names[i][units_from_ui]] = 12 # TODO: add conv factor
                                else:
                                    # TODO: convert
                                    print('convert')

                    entered_list.append((prop_names[i],units_from_ui[i]))




            entered_list = []
            print(snapshot[2])
            m = [i.strip() for i in snapshot[2]]
            snapshot[2] = m
            print(final_md['2_method'])
            print('\nmethod\n')
            for method in final_md['2_method']:
                if method not in snapshot[2] and method not in entered_list:
                    cur.execute('INSERT INTO Model(Method_name) VALUES(%s)',[method])
                    entered_list.append(method)
            
            entered_list = []
            for func in final_md['2_functional']:
                if func not in snapshot[3] and func not in entered_list:
                    cur.execute('INSERT INTO Functional(name) VALUE(%s)',[func])
                    entered_list.append(func)

            entered_list = []
            for basis in final_md['2_basis']:
                if basis not in snapshot[4] and basis not in entered_list:
                    cur.execute('INSERT INTO Basis_set(name) VALUE(%s)',[basis])
                    entered_list.append(basis)

            entered_list = []
            for ff in final_md['2_forcefield']:
                if ff not in snapshot[5] and ff not in entered_list:
                    cur.execute('INSERT INTO Forcefield(name) VALUE(%s)',[ff])
                    entered_list.append(ff)
        else:
            #units = [u[:u.index('(')-1]]
            print('snapshot is false') # DB is empty; executed during first insert
            entered_list = []
            prop_names = [i.lower() for i in prop_names]
            for i in range(len(prop_names)):
                if (prop_names[i],units_from_ui[i]) not in entered_list:
                    if prop_names[i].lower() in unit_list.keys() or prop_type[i] in ['sol','energy']: # If we have this property in our list
                        entered_list.append((prop_names[i],units_from_ui[i]))
                        if 'default' in units_from_ui[i]: # If the unit entered by the user is the default unit for this property in our list
                            u = units_from_ui[i]
                            cur.execute('INSERT INTO Property(Property_str,Unit) VALUES(%s,%s)',[prop_names[i],u[:u.index('(')-1]])
                        else: # if the unit entered by the user is not the default unit in our list
                            if prop_type[i] == 'sol':
                                all_units = list(unit_list['solubility parameters'].keys())
                                def_unit = all_units[0]
                            elif prop_type[i] == 'energy':
                                all_units = list(unit_list['energy'].keys())
                                def_unit = all_units[0]
                            else:
                                all_units = list(unit_list[prop_names[i]].keys())
                                def_unit = all_units[0] # fetching the default unit
                            #print(list(unit_list[prop_names[i]].keys()))
                            if units_from_ui[i] not in all_units: # if we do not have this unit listed for the given property in our list
                                conv_factor = float(units_from_ui[i].split(',')[1])
                                units_from_ui[i] = units_from_ui[i].split(',')[0]
                                if prop_type[i] == 'sol':
                                    unit_list['solubility parameters'][units_from_ui[i]] = conv_factor
                                elif prop_type[i] == 'energy':
                                    unit_list['energy'][units_from_ui[i]] = conv_factor
                                else:
                                    unit_list[prop_names[i]][units_from_ui[i]] = conv_factor
                                unit_flag=True
                            else:
                                if prop_type[i] == 'sol':
                                    conv_factor = unit_list['solubility parameters'][units_from_ui[i]] 
                                elif prop_type[i] == 'energy':
                                    conv_factor = unit_list['energy'][units_from_ui[i]]
                                else:
                                    conv_factor = unit_list[prop_names[i]][units_from_ui[i]]

                            data[props[i]] = data[props[i]] * conv_factor
                            cur.execute('INSERT INTO Property(Property_str,Unit) VALUES(%s,%s)',[prop_names[i],def_unit])
                            print(def_unit)
                            print(unit_list)
                    else: # If the property is not in our list and it is not of the types energy or solubility parameters
                        print(prop_names[i])
                        conv_factor = float(units_from_ui[i].split(',')[1])
                        units_from_ui[i] = units_from_ui[i].split(',')[0]
                        unit_flag = True
                        d = {}
                        d[units_from_ui[i]+' (default)'] = 1.0
                        unit_list[prop_names[i]] = d
                        cur.execute('INSERT INTO Property(Property_str,Unit) VALUES(%s,%s)',[prop_names[i],units_from_ui[i]])
                        # TODO: add property and unit to unit_list and insert into DB
                    entered_list.append((prop_names[i],units_from_ui[i]))

            entered_list = []
            for method in final_md['2_method']:
                if method not in entered_list:
                    cur.execute('INSERT INTO Model(Method_name) VALUES(%s)',[method])
                    entered_list.append(method)
            
            entered_list = []
            for func in final_md['2_functional']:
                if func not in entered_list:
                    cur.execute('INSERT INTO Functional(name) VALUE(%s)',[func])
                    entered_list.append(func)

            entered_list = []
            for basis in final_md['2_basis']:
                if basis not in entered_list:
                    cur.execute('INSERT INTO Basis_set(name) VALUE(%s)',[basis])
                    entered_list.append(basis)

            entered_list = []
            for ff in final_md['2_forcefield']:
                if ff not in entered_list:
                    cur.execute('INSERT INTO Forcefield(name) VALUE(%s)',[ff])
                    entered_list.append(ff)

        
        ## populating the functional table
        #cur.execute("SELECT name FROM Functional")
        #functionals = cur.fetchall()
        #entered_list = []
        #for func in meta_data['2_functional']:
        #    if any(func in i for i in functionals) or func in entered_list:
        #        pass
        #    else:
        #        entered_list.append(func)
        #        cur.execute("INSERT INTO Functional(name) VALUES(%s)",[func])

        #print('entered functional')

        ## populating the basis set table
        #cur.execute("SELECT name FROM Basis_set")
        #basis_sets = cur.fetchall()
        #entered_list = []
        #for basis in meta_data['2_basis']:
        #    if any(basis in i for i in basis_sets) or basis in entered_list:
        #        pass
        #    else:
        #        entered_list.append(basis)
        #        cur.execute("INSERT INTO Basis_set(name) VALUES(%s)",[basis])

        #print('entered basis_set')

        ## populating the forcefield table
        #cur.execute("SELECT name FROM Forcefield")
        #forcefields = cur.fetchall()
        #entered_list = []
        #for ff in meta_data['2_forcefield']:
        #    if any(ff in i for i in forcefields) or ff in entered_list or ff=='None':
        #        pass
        #    else:
        #        entered_list.append(basis)
        #        cur.execute("INSERT INTO Forcefield(name) VALUES(%s)",[ff])
        
        #print('entered forcefield')

        # populating the molecule table

        cur.execute("SELECT "+ ''.join(i.lower()+',' for i in mol_ids.keys())+"MW from Molecule")
        molecules = cur.fetchall()
        new_entries=[]

        # check last id for molecule, add molecule index, melt dataframe, add property and method index using lambda and dictionary
        # check for molecule identifier (name, IUPAC, )
        pybel_identifiers = {'smiles':'smiles','standard_inchi_key':'inchikey','standard_inchi':'inchi'}
        #print(mol_ids)
        
        for mol in range(len(data)):
            mw_flag = False
            row = []
            for k,v in mol_ids.items():
                if k.lower() in pybel_identifiers.keys():
                    try:
                        m = pybel.readstring(pybel_identifiers[k.lower()],data.loc[mol][v])
                        if k.lower() == 'smiles':
                            iden = m.write('can').strip()
                            if mw_flag == False:
                                mw = m.molwt
                                mw = round(mw,3)
                                mw_flag = True
                        else:
                            iden = m.write(pybel_identifiers[k.lower()]).strip()
                            if mw_flag == False:
                                mw = m.molwt
                                mw = round(mw,3)
                                mw_flag = True
                    except:
                        db = db.replace('_',' ')
                        return render_template('temp_insert.html',title=db,all_dbs = all_dbs, err_msg='Invalid SMILES on row number {}.'.format(str(mol)))
                    row.append(iden)
                else:
                    row.append(data.loc[mol][v])
                    if mw_flag == False:
                        mw_flag = True
                        url = 'http://cactus.nci.nih.gov/chemical/structure/'
                        try:
                            url = url + data.loc[mol][v] + '/mw'
                            ans = urlopen(url).read().decode('utf8')
                        except HTTPError:
                            mw = NULL
            if mw_flag == True:
                row.append(mw)
            new_entries.append(row)

        ## temporary fix to enable for case-insensitivity in molecule-identifier
        molecules = [list(x) for x in molecules]
        new_entries = [list(x) for x in new_entries]
        #print(molecules)
        #print(new_entries)
        #for i in molecules:
        #    i[1] = i[1].lower()
        #for i in new_entries:
        #    i[1] = i[1].lower()
        molecules = tuple(tuple(x) for x in molecules)
        new_entries = tuple(tuple(x) for x in new_entries)
        required_entries = list(set(new_entries) - set(molecules))
        if len(mol_ids.keys())>1:
            mol_q = ''.join(i+',' for i in mol_ids.keys())
            vals = ''.join('%s,' for i in range(len(mol_ids)+1))
        else:
            mol_q = list(mol_ids.keys())[0] + ','
            vals = '%s,%s,'
        cur.executemany('INSERT INTO Molecule('+mol_q+'MW) VALUE('+vals[:-1]+')',required_entries)
        print('mol done')
        ### populating the credit table
        ### todo: figure out how to deal with the credit/publication
        ### cur.execute('INSERT INTO %s.Credit(DOI) VALUES(%s)'%db,' ')
        #print(list(mol_ids.values()))
        cols=[i for i in props]
        for i in list(mol_ids.values()):
            cols.append(i)

        data = data[cols]
        ##if smi_col!='':
        ##    cols.append(smi_col)
        ##if mol_identifier!='':
        ##    cols.append(mol_identifier)
        ##data = data[cols]
        # populating the values table

        cur.execute('SELECT id,Property_str from Property')
        all_props = cur.fetchall()
        prop_id = dict(map(reversed,all_props))
        cur.execute('SELECT id,Method_name from Model')
        all_models = cur.fetchall()
        model_id = dict(map(reversed,all_models))
        cur.execute('SELECT id,name from Functional')
        all_functionals = cur.fetchall()
        functional_id=dict(map(reversed,all_functionals))
        cur.execute('SELECT id,name from  Basis_set')
        all_basis = cur.fetchall()
        basis_id = dict(map(reversed,all_basis))
        cur.execute('SELECT id,name from Forcefield')
        all_ff = cur.fetchall()
        ff_id = dict(map(reversed,all_ff))
        mol_q = 'ID,'+ list(mol_ids.keys())[0]
        cur.execute("SELECT "+mol_q+" from Molecule")
        all_mols = cur.fetchall()
        molecule_id = dict(map(reversed,all_mols))
        #print(data.columns)
        #print(molecule_id)
        data['molecule_id']=data[list(mol_ids.values())[0]].apply(lambda a: molecule_id[pybel.readstring('smi',a).write('can').strip()])
        #print(data)
        
        data.drop(mol_ids.values(),1,inplace=True)
        data = data.melt('molecule_id')
        #print(prop_id)
        #print(data)

        data['property_id']=data['variable'].apply(lambda a: prop_id[prop_names[props.index(a)]]) #unable to recognize second column with same property name because it is using lambda to assign property ID
        #print(data.head())
        #print(props)
        print(model_id[final_md['2_method'][props.index(props[0])]])
        data['model_id']=data['variable'].apply(lambda a: model_id[final_md['2_method'][props.index(a)]])
        print(data)
        data['functional_id']=data['variable'].apply(lambda a: functional_id[final_md['2_functional'][props.index(a)]])
        data['basis_id']=data['variable'].apply(lambda a: basis_id[final_md['2_basis'][props.index(a)]])
        data['ff_id']=data['variable'].apply(lambda a: ff_id[final_md['2_forcefield'][props.index(a)]])
        data.drop('variable',1,inplace=True)
        to_drop = []
        id = tuple(data['molecule_id'])
        #print(data)
        cur.execute('select * from Value where molecule_id in {}'.format(str(id)))
        vals = cur.fetchall()
        vals = [list(x) for x in vals]
        vals = pd.DataFrame(vals, columns=['id','value','model_id', 'property_id','molecule_id', 'functional_id','basis_id','ff_id'])
        check_data = data.drop('value',1)
        vals = vals[check_data.columns]
        vals = [list(vals.loc[x]) for x in range(len(vals))]
        for i in range(len(data)):
            if list(check_data.loc[i]) in vals:
                to_drop.append(i)
        data.drop(to_drop,0,inplace=True)
        #print(data)
        if len(data) == 0:
            return render_template('temp_insert.html',title=db,all_dbs=all_dbs,err_msg='Duplicate entries for all molecules exist.')
        else:
            cur.executemany('INSERT INTO Value(molecule_id,num_value,property_id,model_id,functional_id,basis_id,forcefield_id) VALUES(%s,%s,%s,%s,%s,%s,%s)',data.values.tolist())
        db = db.replace('_chembddb','')
        db = db.replace('_',' ')
        if 'download-submit' in final_md:
            df = pd.DataFrame()
            df['Properties'] = [i[0] for i in props]
            df['Units'] = final_md['2_unit']
            df['Method'] = final_md['2_method']
            df['Functional'] = final_md['2_functional']
            df['Basis_set'] = final_md['2_basis']
            df['Forcefield'] = final_md['2_forcefield']
            df['conf'] = df[df.columns].apply(lambda x: ','.join(x), axis = 1)
            confs = ','.join(r for r in df['conf'])
            print(confs)
            cur.execute('INSERT INTO Configuration(ID, conf) VALUES(%s,%s)',[str(0),confs])
        if unit_flag == 'True':
            cur.execute('USE unit_list_chembddb;')
            cur.execute('INSERT INTO Main(id,unit_str) VALUE(1,%s);',[unit_list])
        con.commit()
        if to_drop!=[]:
            all_dbs.append(db)
            return render_template('temp_insert.html',title=db,all_dbs=all_dbs,init='True',err_msg='A few molecules were not entered due to duplicate entries',success_msg='The database has been successfully populated')
        else:
            return render_template('temp_insert.html',title=db,all_dbs=all_dbs,init='True',success_msg='The database has been successfully populated')

    else:
        # default landing page
        return render_template('temp_insert.html',all_dbs=all_dbs,init='True',snapshot='')

@app.route('/search',methods=['GET','POST'])
def search():
    #cur.execute('show databases;') 
    #all_dbs_tup=cur.fetchall()
    #all_dbs=[]
    #for i in all_dbs_tup:
    #    if '_chembddb' in i[0]:
    #        m=i[0]
    #        all_dbs.append((m[:-9],))
    global all_dbs
    return render_template('search.html',all_dbs=all_dbs)

@app.route('/search_db<db>',methods=['GET','POST'])
def search_db(db):
    #global cur,con
    all_dbs=[]
    cur.execute('show databases;')
    all_dbs_tup=cur.fetchall()
    #print(all_dbs_tup)
    for i in all_dbs_tup:
        if '_chembddb' in i[0] and 'unit_list' not in i[0]:
            m=i[0]
            all_dbs.append((m[:-9],))
    db=db[1:-1]
    db=db+'_chembddb'
    #cur.execute('USE INFORMATION_SCHEMA')
    #result=cur.execute('SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA WHERE SCHEMA_NAME=\'%s\''%db)
    #if result == 0:
        #error_message='Database does not exist'
        #return render_template('Search.html',all_dbs=all_dbs,err_msg=error_message)
    cur.execute('USE %s'%db)
    db = db.replace('_chembddb','')
    cur.execute('Select * from Property')
    properties=cur.fetchall()
    #print(properties)
    cur.execute('Select id,method_name from Model')
    results=cur.fetchall()
    cur.execute("Select * from Functional")
    functionals=cur.fetchall()
    cur.execute("Select * from basis_set")
    basis_sets=cur.fetchall()
    cur.execute("Select * from forcefield")
    forcefields=cur.fetchall()
    methods=[]
    # search_results useful only in case of smiles_search because the same results are used ... single db call
    global search_results
    global sql
    global ini
    global fin
    global n_res
    global noprev
    global nonext
    global keys
    global multiprop
    multiprop = False
    for i in results:
        methods.append(i[1])
    if request.method == 'POST' and 'search-query' in request.form:
        # This section is to create the query based on the options provided by the user
        from_form = request.form
        from_form = from_form.to_dict(flat=False)
        keys=[i for i in from_form if '_id' in i]
        min_max_err=False
        min_max_prop=[]
        props=[]
        p = []
        # tuple of tuples to list of tuples
        for pr in properties:
            p.append(list(pr))
        properties = p
        if len(keys)>0:
            sql='select value.molecule_id,molecule.SMILES,model.method_name,functional.name,basis_set.name,forcefield.name,Property.Property_str, value.num_value from molecule inner join Value on molecule.id=value.molecule_id inner join property on property.id=value.property_id inner join model on model.id=value.model_id inner join functional on functional.id=value.functional_id inner join basis_set on basis_set.id = value.basis_id inner join forcefield on forcefield.id=value.forcefield_id where '
            print(from_form)
            for k in keys:
                prop_id=int(from_form[k][0])
                props.append(prop_id)
                from_val=float(from_form[k[:-3]+'_from_val'][0])
                to_val=float(from_form[k[:-3]+'_to_val'][0])
                print(prop_id)
                properties[prop_id-1].append(from_val)
                properties[prop_id-1].append(to_val)
                if from_val > to_val:
                    min_max_err=True
                sql = sql[:sql.rfind('where')+6] + 'molecule_id in (select molecule_id from value where value.property_id={0} and value.num_value>{1} and value.num_value<{2}) and '.format(prop_id,from_val,to_val) + sql[sql.rfind('where')+6:]
            if len(keys)!=0:
                sql=sql[:-5]
            valsid=' and value.property_id in '
            for i in range(len(props)):
                if i > 0:
                    valsid = valsid + ',' + str(props[i])
                else:
                    valsid = valsid + '(' + str(props[i])
            valsid = valsid +')'
            sql = sql +valsid
        else:
            sql = 'select id, SMILES, MW from Molecule where '
        MW_to = None
        if 'MW' in from_form:
            from_val=float(from_form['MW_from_val'][0])
            to_val=float(from_form['MW_to_val'][0])
            MW_from = from_val
            MW_to = to_val
            if from_val > to_val:
                min_max_err=True
            if len(keys)!=0:
                sql=sql+" and molecule.MW > {} and molecule.MW < {} ".format(float(from_form['MW_from_val'][0]),float(from_form['MW_to_val'][0]))
            else:
                sql="select id,SMILES,MW from Molecule where Molecule.MW > {} and Molecule.MW < {} ".format(float(from_form['MW_from_val'][0]),float(from_form['MW_to_val'][0]))
                keys.append('MW')
        if 'smiles_search' in from_form:
            if len(keys)==0:
                sql = 'select id, SMILES from Molecule'

        if 'method' in from_form:
            met_id=0
            for m in results:
                if m[1]==from_form['method_name'][0]:
                    met_id=m[0]
            if len(keys)==0:
                sql=sql+" Value.model_id = {}".format(met_id)
            else:
                sql=sql+" and Value.model_id ={}".format(met_id)

        if 'func' in from_form:
            if len(keys)==0:
                sql=sql+' Value.functional_id={}'.format(from_form['functional_name'][0])
            else:
                sql=sql+' and Value.functional_id={}'.format(from_form['functional_name'][0])
        
        if 'basis' in from_form:
            if len(keys)==0:
                sql=sql+' Value.basis_id={}'.format(from_form['basis_set'][0])
            else:
                sql=sql+' and Value.basis_id={}'.format(from_form['basis_set'][0])
        
        if 'ff' in from_form:
            if len(keys)==0:
                sql=sql+' Value.forcefield_id={}'.format(from_form['forcefield'][0])
            else:
                sql=sql+' and Value.forcefield_id={}'.format(from_form['forcefield'][0])
        
        global counts
        counts_q = 'select count(*) '+sql[sql.find('from'):] +';'

        # because smiles_search will get all results from the db because substructure matching is required
        if ('property' not in sql and 'MW' not in sql) or len(keys)>1:
            counts = -1
        else:
            sql = sql + 'limit 50'
            cur.execute(counts_q)
            counts = cur.fetchall()
            counts = counts[0][0]

        sql=sql+';'
        # query creation ends here
        temp_col=[]              
        temp_met=[]
        if counts == 0:
            if min_max_err==True:
                n_res = 'Min value entered is > Max value entered for one of the fields above.' 
                columns=''
            else:               
                n_res = 'Number of results='+ str(counts)+'. No such candidates exist in your database'
                columns=''
            if 'MW' in from_form:
                return render_template('search_db.html',MW_from=MW_from, MW_to=MW_to,properties=properties,columns=columns,methods=methods,n_res=n_res,basis=basis_sets,functionals=functionals,forcefields=forcefields,all_dbs=all_dbs,title=db)
            else:
                return render_template('search_db.html',properties=properties,columns=columns,methods=methods,n_res=n_res,basis=basis_sets,functionals=functionals,forcefields=forcefields,all_dbs=all_dbs,title=db)
        else:
            # executing the query
            cur.execute(sql)
            data1=cur.fetchall()
            data, columns = post_process(sql, data1)
            # the following section is just to format the column headers that appear on the html page
            for c in columns:
                if '-' in c:
                    temp_col.append(c.split('-')[0])
                    if len(c.split('-'))>2:
                        temp_met.append(c.split('-')[1]+'-'+c.split('-')[2])
                    else:
                        temp_met.append(c.split('-')[1])
                else:
                    if 'MW' in c:
                        temp_met.append('pybel')
                    else:
                        temp_met.append('')
                    temp_col.append(c)
            # substructure matching using pybel
            try:
                smi_val = None
                if 'smiles_search' in from_form:
                    smarts = pybel.Smarts(from_form['smiles'][0])
                    smi_val = smarts
                    for i in range(len(data)):
                        mol = pybel.readstring("smi",data.loc[i]['SMILES'])
                        smarts.obsmarts.Match(mol.OBMol)
                        if len(smarts.findall(mol))==0:
                            data.drop(i,0,inplace=True)
                    if len(data)==0:
                        n_res='Number of results='+ str(len(data))+'\nNo such candidates exist in your database'
                    else:
                        n_res=len(data)
                        counts = n_res
                    if n_res>50:
                        search_results = data
                        search_results.columns = data.columns
                        data = data[:51]
                elif len(keys)>1:
                    if len(data)==0:
                        n_res='Number of results='+ str(len(data))+'\nNo such candidates exist in your database'
                    else:
                        n_res=len(data)
                        counts = n_res
                    if n_res>50:
                        search_results = data
                        search_results.columns = data.columns
                        data = data[:51]
                else:
                    n_res = counts
            except:
                n_res='Invalid Smarts entered'
                data=pd.DataFrame()
            
            desc=['','']
            columns =[]
            # creating tuple of tuples for column headers (required for html page)
            for i in range(len(temp_met)):
                columns.append((temp_col[i],temp_met[i]))
            # calculating statistics for each page
            for i in data.columns[2:]:
                desc.append('mean={}, std={}, min={}, max={}'.format(data[i].describe()['mean'].round(2),data[i].describe()['std'].round(2),data[i].describe()['min'].round(2),data[i].describe()['max'].round(2)))
            data = tuple(data.itertuples(index=False,name=None))
            if len(columns) == 2:
                to_order = False
            else:
                to_order = True
            noprev=True
            db = db.replace('_chembddb','')
            ini=0
            if type(n_res) != str and n_res < 50:
                fin = n_res
                nonext=True
            else:
                nonext=False
                fin = 50
            if 'MW' in from_form:
                return render_template('search_db.html',ini=ini,fin=fin, MW_from=MW_from, MW_to=MW_to,data = data,properties=properties,columns=columns,temp_met=temp_met,methods=methods,n_res=n_res,basis=basis_sets,functionals=functionals,to_order = to_order,forcefields=forcefields,all_dbs=all_dbs,title=db,desc=desc,noprev=noprev,nonext=nonext)
            else:
                return render_template('search_db.html',ini=0,fin=fin,data = data,properties=properties,columns=columns,temp_met=temp_met,methods=methods,n_res=n_res,basis=basis_sets,functionals=functionals,forcefields=forcefields,to_order=to_order, all_dbs=all_dbs,title=db,desc=desc,noprev=noprev,nonext=nonext)
    elif 'next-50' in request.form:
        nonext=False
        from_form = request.form
        from_form = from_form.to_dict(flat=False)
        n_res_done = 0
        if ' offset' in sql:
            # checking the offset in the previous sql query, n_res_done tells us how many results have been displayed already
            n_res_done = int(sql[sql.rfind('offset')+7:-1])
            sql = sql[:sql.rfind('offset')+6] + ' '+str(n_res_done + 50) +';'
            n_res_done = n_res_done + 50
        else:
            n_res_done = 50
            sql = sql[:-1]+' offset '+str(n_res_done)+';'
        if ('property' not in sql and 'MW' not in sql) or len(keys)>1:
            data = search_results[n_res_done:n_res_done+51]
            columns = list(data.columns)
            if len(keys)==0:
                to_order = False
            else:
                columns=[c.replace('(NA/NA)','') for c in columns]
                columns=[c.replace('(na/na)','') for c in columns]
                columns=[c.replace('(NA)','') for c in columns]
                columns=[c.replace('(na)','') for c in columns]
                to_order = True
        else:   
            to_order = True         
            cur.execute(sql)
            data1=cur.fetchall() 
            data, columns = post_process(sql, data1)
        temp_col=[]
        temp_met=[]
        for c in columns:
            if '-' in c:
                temp_col.append(c.split('-')[0])
                if len(c.split('-'))>2:
                    temp_met.append(c.split('-')[1]+'-'+c.split('-')[2])
                else:
                    temp_met.append(c.split('-')[1])
            else:
                temp_col.append(c)
                if 'MW' in c:
                    temp_met.append('pybel')
                else:
                    temp_met.append('')

        desc=['','']
        columns =[]
        ini = n_res_done
        if (counts - n_res_done) < 50: 
            fin = counts
            nonext = True
        else:
            nonext = False
            fin = n_res_done + 50                
        noprev=False
        if temp_met!=[]:
            for i in range(len(temp_met)):
                columns.append((temp_col[i],temp_met[i]))
            for i in data.columns[2:]:
                desc.append('mean={}, std={}, min={}, max={}'.format(data[i].describe()['mean'].round(2),data[i].describe()['std'].round(2),data[i].describe()['min'].round(2),data[i].describe()['max'].round(2)))
        if 'order by' in sql:
            if 'DESC' in sql:        
                data=data.sort_values(by=data.columns[-1],ascending=False)
            else:
                data=data.sort_values(by=data.columns[-1])
        data = tuple(data.itertuples(index=False,name=None))
        if 'MW' in sql and 'value.property_id=' not in sql:
            return render_template('search_db.html',ini=ini,fin=fin, to_order=to_order,noprev=False,nonext=nonext,data = data,properties=properties,columns=columns,temp_met=temp_met,methods=methods,n_res=n_res,basis=basis_sets,functionals=functionals,forcefields=forcefields,all_dbs=all_dbs,title=db,desc=desc)
        else:
            return render_template('search_db.html',ini=ini,fin=fin,to_order=to_order,nonext=nonext,noprev=False,data = data,properties=properties,columns=columns,temp_met=temp_met,methods=methods,n_res=n_res,basis=basis_sets,functionals=functionals,forcefields=forcefields,all_dbs=all_dbs,title=db,desc=desc)
    elif 'prev-50' in request.form:
        noprev=False
        from_form = request.form
        from_form = from_form.to_dict(flat=False)
        n_res_done = 0
        if ' offset' in sql:
            # checking the offset in the previous sql query, this number tells us how many results have been displayed already
            n_res_done = int(sql[sql.rfind('offset')+7:-1]) - 50
            sql = sql[:sql.rfind('offset')+6] + ' '+str(n_res_done) +';'
            noprev = False
        if n_res_done ==0:
            noprev = True
            nonext = False
        if ('property' not in sql and 'MW' not in sql) or len(keys)>1:
            data = search_results[n_res_done:n_res_done+51]
            columns = list(data.columns)
            if len(keys)==0:
                to_order = False
            else:
                columns=[c.replace('(NA/NA)','') for c in columns]
                columns=[c.replace('(na/na)','') for c in columns]
                columns=[c.replace('(NA)','') for c in columns]
                columns=[c.replace('(na)','') for c in columns]
                to_order = True
        else:
            cur.execute(sql)
            data1=cur.fetchall()
            data, columns = post_process(sql, data1)
            to_order = True
        temp_col=[]
        temp_met=[]

        for c in columns:
            if '-' in c:
                temp_col.append(c.split('-')[0])
                if len(c.split('-'))>2:
                    temp_met.append(c.split('-')[1]+'-'+c.split('-')[2])
                else:
                    temp_met.append(c.split('-')[1])
            else:
                temp_col.append(c)
                if 'MW' in c:
                    temp_met.append('pybel')
                else:
                    temp_met.append('')

        desc=['','']
        columns =[]
        fin = n_res_done + 50     
        if temp_met!=[]:
            for i in range(len(temp_met)):
                columns.append((temp_col[i],temp_met[i]))
            for i in data.columns[2:]:
                desc.append('mean={}, std={}, min={}, max={}'.format(data[i].describe()['mean'].round(2),data[i].describe()['std'].round(2),data[i].describe()['min'].round(2),data[i].describe()['max'].round(2)))
        
        if 'order by' in sql:
            if 'DESC' in sql:        
                data=data.sort_values(by=data.columns[-1],ascending=False)
            else:
                data=data.sort_values(by=data.columns[-1])
        data = tuple(data.itertuples(index=False,name=None))
        ini = n_res_done
        if 'MW' in sql and 'value.property_id=' not in sql:
            return render_template('search_db.html',ini=n_res_done,fin=fin,to_order=to_order, noprev = noprev, nonext=False,data = data,properties=properties,columns=columns,temp_met=temp_met,methods=methods,n_res=n_res,basis=basis_sets,functionals=functionals,forcefields=forcefields,all_dbs=all_dbs,title=db,desc=desc)
        else:
            return render_template('search_db.html',ini=n_res_done,fin=fin, to_order=to_order, noprev = noprev,nonext=False,data = data,properties=properties,columns=columns,temp_met=temp_met,methods=methods,n_res=n_res,basis=basis_sets,functionals=functionals,forcefields=forcefields,all_dbs=all_dbs,title=db,desc=desc)        
    elif 'download_csv' in request.form or 'download_json' in request.form:
        from_form = request.form
        from_form = from_form.to_dict(flat=False)
        desc=['','']
        ## re-executing query to get all results

        if 'MW' not in sql and 'property' not in sql:
            data = search_results
            columns = search_results.columns
            to_order = False
        else:
            if 'MW' in sql and 'property' not in sql:
                sql = sql[:sql.rindex('limit')]+';'
            else:
                sql = sql[:sql.rindex(')')+1]+';'                   
            cur.execute(sql)        
            all_results = cur.fetchall()
            data, columns = post_process(sql, all_results)
            to_order = True
        
        if 'download_json' in request.form:
            import json
            data.to_json('results.json')
            msg='Results have been downloaded as results.json'
        else:
            print(data)
            data.to_csv('results.csv',index=None)
            msg='Results have been downloaded as results.csv'
        ## results that go to the html are still limited to 50
        for i in data.columns[2:]:
            desc.append('mean={}, std={}, min={}, max={}'.format(data[i].describe()['mean'].round(2),data[i].describe()['std'].round(2),data[i].describe()['min'].round(2),data[i].describe()['max'].round(2)))
        # search_results.to_csv('results.csv',index=None)
        
        columns=[]
        for i in data.columns:
            if '-' not in i:
                if 'MW' in i:
                    columns.append((i,'pybel'))
                else:
                    columns.append((i,''))
            else:
                if len(i.split('-')) > 2:
                    columns.append((i.split('-')[0],i.split('-')[1]+'-'+i.split('-')[2]))
                else:
                    columns.append((i.split('-')[0],i.split('-')[1])) 
        data = tuple(data.itertuples(index=False,name=None))
        return render_template('search_db.html',data = data,ini=ini, fin=fin, to_order=to_order, properties=properties,columns=columns,methods=methods,msg=msg,n_res=n_res,functionals=functionals,basis=basis_sets,forcefields=forcefields,all_dbs=all_dbs,noprev=noprev,nonext=nonext,title=db,desc=desc)
    elif 'orderby_property' in request.form:
        from_form=request.form
        from_form = from_form.to_dict(flat=False)
        ascending = True
        if 'ascending' in from_form['select_order']:
            if 'order by' not in sql:
                if 'property' not in sql and 'MW' in sql:
                    sql = sql[:sql.rindex('limit')] + ' order by molecule.MW ' +sql[sql.rindex('limit'):]
                else:
                    sql = sql[:sql.rindex(')')+1]+ ' order by Value.num_value ' + sql[sql.rindex(')')+1:]
            elif 'order by' in sql and 'DESC' in sql:
                if 'property' not in sql and 'MW' in sql:
                    sql = sql[:sql.rindex('DESC')] + sql[sql.rindex('DESC')+5:]
                else:
                    sql = sql[:sql.rindex('value')+5] + ' ' + sql[sql.rindex('value')+11:]
        else:
            ascending = False
            if 'order by' in sql and 'DESC' not in sql:
                if 'property' not in sql and 'MW' in sql:
                    sql = sql[:sql.rindex('MW')+3] + 'DESC ' + sql[sql.rindex('MW')+3:]   
                else:               
                    sql = sql[:sql.rindex('value')+5] + ' DESC' + sql[sql.rindex('value')+5:]
            elif 'order by' not in sql:
                if 'property' not in sql and 'MW' in sql:
                    sql = sql[:sql.rindex('limit')] + ' order by molecule.MW DESC ' +sql[sql.rindex('limit'):]
                else:
                    sql = sql[:sql.rindex(')')+1]+ ' order by value.num_value DESC ' + sql[sql.rindex(')')+1:]
        
        if sql.count('value.num_value') > 4:
            multiprop = True
            sql = sql[:sql.rfind(')')+1]+';'
        cur.execute(sql)
        all_results = cur.fetchall()
        data, columns = post_process(sql, all_results)
        search_results = data              
        search_results.columns = columns
        if 'MW' not in sql and 'property' in sql:
            search_results = search_results.sort_values(by=from_form['property_orderby'], ascending = ascending)
        desc=['','']
        data = tuple(search_results[:50].itertuples(index=False,name=None))
        for i in search_results.columns[2:]:
            desc.append('mean={}, std={}, min={}, max={}'.format(search_results[i].describe()['mean'].round(2),search_results[i].describe()['std'].round(2),search_results[i].describe()['min'].round(2),search_results[i].describe()['max'].round(2)))
        columns=[]
        for i in search_results.columns:
            if '-' not in i:
                if 'MW' in i:
                    columns.append((i,'pybel'))
                else:
                    columns.append((i,''))
            else:
                if len(i.split('-')) > 2:
                    columns.append((i.split('-')[0],i.split('-')[1]+'-'+i.split('-')[2]))
                else:
                    columns.append((i.split('-')[0],i.split('-')[1]))
        return render_template('search_db.html',data = data,properties=properties,to_order = True, ini=ini,fin=fin,noprev=noprev,nonext=nonext,columns=columns,methods=methods,n_res=n_res,basis=basis_sets,functionals=functionals,forcefields=forcefields,all_dbs=all_dbs,title=db,desc=desc)

    else:
        return render_template('search_db.html',properties=properties,methods=methods,basis=basis_sets, functionals=functionals, forcefields=forcefields,all_dbs=all_dbs,title=db)

@app.route('/molecule-<dbid>',methods=['GET','POST'])
def molecule(dbid):
    # cur,conn=connect_mysql()
    import urllib.parse
    db=dbid.split('-')[0]
    db=db.replace(' ','_')
    id=dbid.split('-')[1]
    sql = 'SELECT Molecule.id, Molecule.MW, Molecule.SMILES,Molecule.Standard_Inchi,Molecule.Standard_Inchi_Key,CAS_Registry_Number,IUPAC_Name,Other_name,Chemical_Formula,Property.Property_str, Property.Unit, Model.method_name,  Functional.name, Basis_set.name,forcefield.name, Value.num_value from Molecule inner join value on Molecule.id=Value.Molecule_id inner join Property on Property.id=VALUE.property_id INNER JOIN Model on Value.model_id=Model.id inner join functional on functional.id=value.functional_id inner join basis_set on basis_set.id=value.basis_id inner join forcefield on forcefield.id=value.forcefield_id where Molecule.id={}'.format(id)
    db=db+'_chembddb'
    cur.execute('USE {};'.format(db))
    cur.execute(sql)
    result=cur.fetchall()
    mol_data=pd.DataFrame(list(result),columns=['ID','MW','SMILES','Standard_Inchi','Standard_Inchi_Key','CAS_Registry_Number','IUPAC_Name','Other_name','Chemical_Formula','Property','Unit','Method','Functional','Basis_set','Forcefield','Value'])
    mol_data['ALL']=mol_data['ID'].astype(str) +',;'+mol_data['MW'].astype(str)+',;'+mol_data['SMILES']+',;'+mol_data['Standard_Inchi']+',;'+mol_data['Standard_Inchi_Key']+',;'+mol_data['CAS_Registry_Number']+',;'+mol_data['IUPAC_Name']+',;'+mol_data['Other_name']+',;'+mol_data['Chemical_Formula']

    mol_data['Property(Unit)']=mol_data['Property']+' ('+mol_data['Unit']+')\n'+'- '+mol_data['Method']+'('+mol_data['Functional']+'/'+mol_data['Basis_set']+')('+mol_data['Forcefield']+')'
    mol_data=mol_data[['ALL','Property(Unit)','Value']]
    mol_data=mol_data.pivot(index='ALL',columns='Property(Unit)')
    mol_data=mol_data['Value'].reset_index()
    mol_data[['ID','MW','SMILES','Standard_Inchi','Standard_Inchi_Key','CAS_Registry_Number','IUPAC_Name','Other_name','Chemical_Formula']]=mol_data['ALL'].str.split(',;',expand=True)
    #print(url_smi)
    #print(mol_data['SMILES'][0])
    # converting smiles to chemspider ID NOTE: use this for other mol identifiers
    
    #url = "http://cactus.nci.nih.gov/chemical/structure/{}/chemspider_id".format(mol_data['SMILES'][0])
    #print(url)

    cols=['ID','MW','SMILES','Standard_Inchi','Standard_Inchi_Key','CAS_Registry_Number','IUPAC_Name','Other_name','Chemical_Formula']
    for i in mol_data.columns[1:-9]:
        cols.append(i)
    mol_data=mol_data[cols]
    msg = ''

    mol_data,added = molidentfiers.populate_molidentifiers(mol_data)
    mol_data = mol_data.replace('NONE',np.nan)
    mol_data.dropna(axis=1,inplace=True)

    url_smi = urllib.parse.quote_plus(mol_data['SMILES'][0])
    smi = str(mol_data.loc[0]['SMILES'])
    mol_ob = pybel.readstring("smi",smi)
    mymol = pybel.readstring("smi", mol_ob.write("can"))
    mymol.make3D(forcefield='mmff94', steps=50)
    mymol.write('xyz',app.config['UPLOAD FOLDER']+'/chembddb_{}.xyz'.format(mol_data['ID'][0]),overwrite=True)
    with open(app.config['UPLOAD FOLDER']+'/chembddb_{}.xyz'.format(mol_data['ID'][0])) as f:
        xyz = f.read()
        mol_data['XYZ'] = xyz

    added = list(set(added).intersection(set(mol_data.columns)))
    button = []
    for i in range(len(mol_data.columns)):
        if mol_data.columns[i] in added:
            button.append('True')
        else:
            button.append('')
    button = pd.Series(button,index=mol_data.columns,name=1)
    mol_data = mol_data.append(button)
    cols = mol_data.columns[:-1]
    if request.method == 'POST' and 'addtodb' in request.form:
        print('here')
        print(request.form['addtodb'])
        identifier = request.form['addtodb'].split(',')[0]
        value = request.form['addtodb'][len(identifier)+1:]
        if 'InChI=' in value:
            value = value[6:]
        print(value)
        print(identifier)
        q = 'Update Molecule set {}="{}" where id={};'.format(identifier,value,id)
        print(q)
        cur.execute(q)
        con.commit()
        print('done')
        mol_data[identifier][1] = ''
    #ids = molidentifiers()
    if request.method == 'POST' and 'Download' in request.form:
        mol_data.to_json('molecule.json')
        msg = 'Downloaded molecule.json'
    cols=[c.replace('(na/na)','') for c in cols]
    cols=[c.replace('(na)','') for c in cols]
    cols=[c.replace('(na)','') for c in cols]
    mol_data = tuple(mol_data.itertuples(index=False,name=None))
    mol_data = (tuple(cols[:]),)+mol_data
    mol_data=tuple(zip(*mol_data))
    db=db.replace('_',' ')

    return render_template('molecule.html',mol_data=mol_data,columns=cols,title=db,all_dbs=all_dbs,url_smi=url_smi,msg=msg)

@app.route('/delete',methods=['GET','POST'])
def delete(host='',user='',pw='',db=''):
    """
    Delete a database that was created using chembddb or delete data from one.

    Parameters
    ----------
    host: str default=''
        the hostname is the domain name or server name
    user: str default=''
        the username for mysql
    pw: str default=''
        the password for mysql
    db: str default=''
        the name of the database that needs to be set up

    """
    global all_dbs
    if db !='':
        try:
            db = db+'_chembddb'
            a,b = connect_mysql(host=host,user=user,pw=pw)
            cur.execute('drop database %s;'%db)
            return 'Successfully deleted the database'
        except:
            return 'Failed! database does not exist'
    else:
        cur.execute('show databases;')
        all_dbs=[]
        all_dbs_tup=cur.fetchall()
        for i in all_dbs_tup:
            if '_chembddb' in i[0] and 'unit_list' not in i[0]:
                m=i[0]
                all_dbs.append((m[:-9],))
        details=request.form
        details=details.to_dict(flat=True)
        if 'dbname' in details:
            dbname=details['dbname']+'_chembddb'
            cur.execute('use {};'.format(dbname))
            cur.execute('Select * from Property')
            properties=cur.fetchall()
            cur.execute('Select id,method_name from Model')
            results=cur.fetchall()
            cur.execute("Select * from Functional")
            functionals=cur.fetchall()
            cur.execute("Select * from Basis_set")
            basis_sets=cur.fetchall()
            cur.execute("Select * from Forcefield")
            forcefields=cur.fetchall()
            methods=[]
            for i in results:
                methods.append(i[1])
            dbname = dbname[:-9]
        if 'submit' in details:
            return render_template('delete.html',data=True,dbname=dbname,properties=properties,methods=methods,functionals=functionals,basis=basis_sets,forcefields=forcefields,all_dbs=all_dbs)
        elif 'search-query' in details:

            #if 'exampleRadios' not in details:
            if details['exampleRadios'] == 'option_null':
                dbname = dbname+'_chembddb'
                cur.execute('drop database {}'.format(dbname))
                cur.execute('show databases;')
                all_dbs_tup=cur.fetchall()
                all_dbs=[]
                for i in all_dbs_tup:
                    if '_chembddb' in i[0] and 'unit_list' not in i[0]:
                        m=i[0]
                        all_dbs.append((m[:-9],))
                #return render_template('delete.html',data=True,properties=properties,methods=methods,functionals=functionals,basis=basis_sets,forcefields=forcefields,all_dbs=all_dbs,success_msg='database {} deleted'.format(details['dbname']))
                return render_template('delete.html',all_dbs=all_dbs, success_msg = 'database {} deleted successfully'.format(dbname.replace('_chembddb','')))
            elif details['exampleRadios']=='option1':
                if 'MW' in details and details['MW']!='':
                    mw_from=details['MW_from_val']
                    mw_to=details['MW_to_val']
                if 'smiles_search' in details:
                    smi = details['smiles']
                    smi_obj = pybel.readstring('smi',smi)
                    can_smi = smi_obj.write('can').strip()
                    mol_wt = round(smi_obj.molwt,3)
                    cur.execute('select id,SMILES,MW from Molecule where SMILES=\'{0}\' and (MW-{1}) < 0.00001;'.format(can_smi,mol_wt))
                    to_delete=list(cur.fetchall())
                    sql = 'delete from Value where molecule_id={};'.format(to_delete[0][0])
                    cur.execute(sql)
                    cur.execute('delete from Molecule where id={};'.format(to_delete[0][0]))
                    con.commit()
            else:
                dbname = details['dbname'] + '_chembddb'
                cur.execute('USE {};'.format(dbname))
                # name of property check boxes have '_id' in them
                keys=[i for i in details if '_id' in i]
                # find molecule id if smiles_search in details
                # find MW from and to
                for k in keys:
                    prop_id=int(details[k][0])
                    # Remove the property from the database
                    if details[k[:-3]+'_from_val']=='' and details[k[:-3]+'_to_val']=='':
                        sql='DELETE FROM Property WHERE id={};'.format(prop_id)
                        cur.execute(sql)
                        con.commit()
                        cur.execute('use {};'.format(dbname))
                        cur.execute('Select * from Property')
                        properties=cur.fetchall()
                        cur.execute('Select id,method_name from Model')
                        results=cur.fetchall()
                        cur.execute("Select * from Functional")
                        functionals=cur.fetchall()
                        cur.execute("Select * from Basis_set")
                        basis_sets=cur.fetchall()
                        cur.execute("Select * from Forcefield")
                        forcefields=cur.fetchall()
                        methods=[]
                        for i in results:
                            methods.append(i[1])

                    else:
                        from_val=float(details[k[:-3]+'_from_val'])
                        to_val=float(details[k[:-3]+'_to_val'])
                        if from_val > to_val:
                            return render_template('delete.html',data=True, dbname=details['dbname'].replace('_chembddb',''),properties=properties,methods=methods,functionals=functionals,basis=basis_sets,forcefields=forcefields,all_dbs=all_dbs,err_msg='Minimum value for one of the properties is greater than the maximum value for it.')
                        else:
                            sql='DELETE FROM Value WHERE property_id={} and num_value > {} and num_value < {};'.format(prop_id,from_val,to_val)
                            cur.execute(sql)
                            #cur.execute('select id from Molecule')
                            #mol_ids = cur.fetchall()
                            #cur.execute('select molecule_id from Value')
                            #mol_ids_val = cur.fetchall()
                            #to_delete = []
                            #mol_ids_val = set(mol_ids_val)
                            #for i in mol_ids:
                            #    if i not in mol_ids_val:
                            #        to_delete.append(i[0])
                            #print('HERE:\n')
                            #print(to_delete)
                            #cur.execute('delete from Molecule where id={};'.format(str(tuple(to_delete))))
                            con.commit()
            return render_template('delete.html',data=True, dbname=details['dbname'].replace('_chembddb',''),properties=properties,methods=methods,functionals=functionals,basis=basis_sets,forcefields=forcefields,all_dbs=all_dbs,success_msg='Deleted from database {}.'.format(details['dbname'].replace('_chembddb','')))
        else:
            return render_template('delete.html',all_dbs=all_dbs)

def backup_restore(host='',user='',pw='',db='',filename=''):
    import subprocesses
    if db !='':
        pass
    else:
        try:
            s= 'mysqldump -h localhost -P 3306 -u '+user+' -p'+pw+ ' '+db+' --single-transaction  > ' + filename
            subprocesses.Popen(s,shell=True)
            return 'done'
        except:
            return 'error'
        
@app.route('/uploads/<filename>')
def uploaded_file(filename):
    return send_from_directory(app.config['UPLOAD FOLDER'],filename)

def run_config():
    print('Open http://localhost:5000/connect')
    app.run(debug=True)
    for i in os.listdir():
        if 'chembddb' in i and 'xyz' in i:
            os.remove(i)
    # os.rmdir('./uploads')