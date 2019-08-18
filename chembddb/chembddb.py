from flask import Flask, render_template, url_for, request,redirect
import pymysql
import os
import sys
import pandas as pd
from copy import deepcopy
import pybel
from flask import send_from_directory
import numpy as np
import time


app = Flask(__name__)
upload_directory=os.getcwd()
app.config['UPLOAD FOLDER']=upload_directory

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
    global cur,all_dbs
    """establishes mysql connection based on the credentials provided during setup; any changes in credentials should be made in credentails.dat

    Parameters
    ----------
    host=
    user=
    pw=
    Returns
    -------
    cur: cursor object
        pointer to the MySQL database; used to execute all SQL queries
    """ 
    if request.method=='POST':
        cred=request.form
        cred = cred.to_dict(flat=False)
        cur,all_dbs = connect_mysql(host = cred['host'][0], user=cred['username'][0], pw = cred['password'][0])
        # con = pymysql.connect(host = cred['host'][0], user=cred['username'][0], password = cred['password'][0])
        # cur = con.cursor()
        if cur == 'invalid' and all_dbs == 'credentials':
            return render_template('connect.html',err_msg='Invalid Credentials. Did not connect to MySQL')
        else:
            return render_template('connect.html',success_msg='Connection Established',host=cred['host'][0],user=cred['username'][0],password=cred['password'][0],all_dbs=all_dbs)
    else:
        return render_template('connect.html')
    
@app.route('/setup',methods=['GET','POST'])
def setup(host='',user='',pw='',db=''):
    """
    Function to setup the database with the chembddb schema

    Parameters
    ----------
    host='': 
    user='': 
    pw='': 
    name='': 

    Returns
    -------
    True/False in case of success/faliure
    """
    if host != '':
        # print('here')
        b, a = connect_mysql(host=host,user=user,pw=pw)
        db = db +'_chembddb'
    elif request.method=='POST':
        db_details=request.form
        db_details=db_details.to_dict(flat=False)
        db=db_details['dbname'][0]+'_chembddb'
    else:
        all_dbs=[]
        cur.execute('show databases;')
        all_dbs_tup=cur.fetchall()
        for i in all_dbs_tup:
            if '_chembddb' in i[0]:
                m=i[0]
                all_dbs.append((m[:-9],))
        return render_template('setup.html',all_dbs=all_dbs)
    
    all_dbs=[]
    cur.execute('show databases;')
    all_dbs_tup=cur.fetchall()
    # print(all_dbs_tup)
    for i in all_dbs_tup:
        if '_chembddb' in i[0]:
            m=i[0]
            all_dbs.append((m[:-9],))
    # print(all_dbs)
    cur.execute('USE INFORMATION_SCHEMA')
    result=cur.execute('SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA WHERE SCHEMA_NAME=\'%s\''%db)
    # result=cur.fetchall()
    # print(result)
    if result == 0:
        cur.execute('CREATE DATABASE %s;'%db)
        cur.execute('USE %s;'%db)
        cur.execute('CREATE TABLE `%s`.`Molecule` (`id` INT NOT NULL AUTO_INCREMENT,`SMILES_str` VARCHAR(500) DEFAULT \'NONE\', `Molecule_identifier` VARCHAR(200) DEFAULT \'NONE\',`MW` FLOAT, PRIMARY KEY (`id`));'%db)
        # cur.execute('CREATE TABLE `%s`.`Credit`(`id` INT NOT NULL AUTO_INCREMENT,`DOI` VARCHAR(100) UNIQUE DEFAULT \'None\',`details` VARCHAR(100) DEFAULT \'None\',PRIMARY KEY (`id`));'%db)
        cur.execute('CREATE TABLE `%s`.`Property`(`id` INT NOT NULL AUTO_INCREMENT,`Property_str` VARCHAR(100) NOT NULL UNIQUE,`Unit` VARCHAR(100) NOT NULL,PRIMARY KEY (`id`));'%db)
        cur.execute('CREATE TABLE `%s`.`Model`(`id` INT NOT NULL AUTO_INCREMENT,`method_name` VARCHAR(100) NOT NULL UNIQUE,`options` VARCHAR(500),PRIMARY KEY (`id`));'%db)
        cur.execute('CREATE TABLE `%s`.`Functional`(`id` INT NOT NULL AUTO_INCREMENT,`name` VARCHAR(100) DEFAULT \'NONE\',PRIMARY KEY (`id`));'%db)
        cur.execute('CREATE TABLE `%s`.`Basis_set`(`id` INT NOT NULL AUTO_INCREMENT,`name` VARCHAR(100) DEFAULT \'NONE\',PRIMARY KEY (`id`));'%db)
        cur.execute('CREATE TABLE `%s`.`Forcefield`(`id` INT NOT NULL AUTO_INCREMENT,`name` VARCHAR(100) DEFAULT \'NONE\',PRIMARY KEY (`id`));'%db)
        # cur.execute('CREATE TABLE `%s`.`Topology`(`id` INT NOT NULL AUTO_INCREMENT,`geometry` VARCHAR(100) NOT NULL,`symbols` VARCHAR(100),`method` VARCHAR(100),`steps` INT,PRIMARY KEY (`id`));'%db)
        cur.execute('CREATE TABLE `%s`.`Value`(`id` INT NOT NULL AUTO_INCREMENT,`num_value` FLOAT NOT NULL,`model_id` INT NOT NULL,`property_id` INT NOT NULL,`molecule_id` INT NOT NULL,`functional_id` INT, `basis_id` INT,`forcefield_id` INT,PRIMARY KEY (`id`));'%db)

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
            if '_chembddb' in i[0]:
                m=i[0]
                all_dbs.append((m[:-9],))
        # cur.ex ecute('ALTER TABLE `%s`.`Value` ADD CONSTRAINT `Value_fk3` FOREIGN KEY (`credit_id`) REFERENCES `Credit`(`id`);'%db)
        if host == '':
            return render_template('setup.html',dbname=db,all_dbs=all_dbs,success_msg='The database has been created.')
        else:
            return True
    else:
        if host == '':
            # print(all_dbs)
            return render_template('setup.html',dbname=db,all_dbs=all_dbs,err_msg='Database already exists.')
        else:
            return False

@app.route('/insert',methods=['GET','POST'])
def insert(host='',user='',pw='',db='',smi_col='',mol_identifier='',conf_file='',data_file=''):
    # db=db[1:-1]
    mi_cols=[]
    # cur,conn=connect_mysql()
    if host !='':
        b,a = connect_mysql(host=host, user=user,pw=pw)
        db=db+'_chembddb'
    elif request.method=='POST':
        config_options=request.form
        # print(config_options)
        config_options=config_options.to_dict(flat=False)
        db=config_options['dbname'][0]
        # print(config_options)
        db = db+'_chembddb'
    else:
        cur.execute('show databases;') 
        all_dbs_tup=cur.fetchall()
        all_dbs=[]
        for i in all_dbs_tup:
            if '_chembddb' in i[0]:
                m=i[0]
                all_dbs.append((m[:-9],))
        return render_template('insert.html',all_dbs=all_dbs)
    
    cur.execute('USE INFORMATION_SCHEMA')
    result=cur.execute('SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA WHERE SCHEMA_NAME=\'%s\''%db)
    if result == 0:
        error_message='Database does not exist'
        if host == '':
            return render_template('insert.html',all_dbs=all_dbs,err_msg=error_message)
        else:
            return False
    else:
        # print('here1')
        # print(files)
        # print(conf_file)

        if type(conf_file) is str and conf_file=='':
            # print('1')
            files=request.files
            files=files.to_dict(flat=False)
            conf_file=files['config_file'][0]
            data_file=files['data_file'][0]
            smi_col=config_options['smiles'][0]
            mol_identifier=config_options['molecule_identifier'][0]
            # db = db+'_chembddb'
            if conf_file.filename=='' or conf_file.filename.rsplit('.',1)[1]!='csv':
                db.replace('_chembddb','')
                db=db.replace('_',' ')
                return render_template('insert.html',title=db,err_msg='No config file provided or incorrect file format. (csv required)')
            elif data_file.filename=='' or data_file.filename.rsplit('.',1)[1]!='csv':
                db.replace('_chembddb','')
                db=db.replace('_',' ')
                return render_template('insert.html',title=db,err_msg='No data file provided or incorrect file format. (csv required)')
            elif smi_col=='' and mol_identifier=='':
                db.replace('_chembddb','')
                db=db.replace('_',' ')
                return render_template('insert.html',title=db,err_msg='No molecule identifiers provided.')
            conf=pd.read_csv(conf_file)
            data=pd.read_csv(data_file)

        elif type(conf_file) is str and conf_file!='':
            # print('2')
            conf=pd.read_csv(conf_file)
            data=pd.read_csv(data_file)
        else:
            # print('3')
            conf = conf_file    
            data = data_file

        conf.replace(np.nan,'na',inplace=True)
        all_prop=True
        all_mols=True
        for prop in conf['properties']:
            if prop not in data.columns:
                all_prop=False
        if all_prop==False:
            db=db.replace('_',' ')
            if host=='':
                return render_template('insert.html',title=db,err_msg='Property(s) in config file do not exist in data.')
            else:
                return False
        else:
            if smi_col!='' and smi_col not in data.columns:
                all_mols=False
            if mol_identifier!='' and mol_identifier not in data.columns:
                all_mols=False
            if all_mols==False:
                if host=='':
                    db=db.replace('_chembddb','')
                    db=db.replace('_',' ')
                    return render_template('insert.html',title=db,err_msg='Identifier(s) listed do not exist in data.')
                else:
                    return False
            else:
                cur.execute('USE %s;'%db)
                db = db.replace('_chembddb','')
                # loop throught he CSV file, check if the smiles value is in the table, if yes, fetch the corresponding id, same goes for property, same goes for method for that property, if it does not exist, fetch the last id and create a new entry

                # populating and property table

                entered_list=[]
                cur.execute("SELECT Property_str from Property")
                properties = cur.fetchall()
                # print(properties)
                for prop, units in zip(conf['properties'],conf['units']):
                    if any(prop in i for i in properties) or prop in entered_list:
                        pass
                    else:
                        entered_list.append(prop)
                        cur.execute("INSERT INTO Property(Property_str,Unit) VALUES(%s,%s)",[prop,units])
                print('property table populated')
                # populating the model table
                cur.execute("SELECT Method_name from Model")
                models = cur.fetchall()
                entered_list=[]
                for method in conf['methods']:
                    if any(method in i for i in models) or method in entered_list:
                        pass
                    else:
                        entered_list.append(method)
                        cur.execute("INSERT INTO Model(Method_name) VALUES(%s)",[method])
                print('method table populated')
                # populating the functional table
                cur.execute("SELECT name FROM functional")
                functionals = cur.fetchall()
                entered_list = []
                for func in conf['functional']:
                    if any(func in i for i in functionals) or func in entered_list:
                        pass
                    else:
                        entered_list.append(func)
                        cur.execute("INSERT INTO Functional(name) VALUES(%s)",[func])
                print('functional table populated')
                # populating the basis_set table
                cur.execute("SELECT name FROM Basis_set")
                basis_sets = cur.fetchall()
                entered_list = []
                for basis in conf['basis']:
                    if any(basis in i for i in basis_sets) or basis in entered_list:
                        pass
                    else:
                        entered_list.append(basis)
                        cur.execute("INSERT INTO Basis_set(name) VALUES(%s)",[basis])

                print('basis table populated')
                # populating the forcefield table
                cur.execute("SELECT name FROM Forcefield")
                forcefields = cur.fetchall()
                entered_list = []
                for ff in conf['forcefield']:
                    if any(ff in i for i in forcefields) or ff in entered_list or ff=='None':
                        pass
                    else:
                        entered_list.append(basis)
                        cur.execute("INSERT INTO Forcefield(name) VALUES(%s)",[ff])
                print('ff table populated')
                # populating the molecule table
                cur.execute("SELECT SMILES_str,Molecule_identifier,MW from Molecule")
                molecules = cur.fetchall()
                new_entries=[]
                # check last id for molecule, add molecule index, melt dataframe, add property and method index using lambda and dictionary
                # check for molecule identifier (name, IUPAC, )
                
                for mol in range(len(data)):
                    if smi_col=='':
                        if data.loc[mol][mol_identifier] not in new_entries:
                            new_entries.append(('None',data.loc[mol][mol_identifier],'None'))
                    elif mol_identifier=='':
                        if data.loc[mol][smi_col] not in new_entries:
                            try:
                                m = pybel.readstring("smiles",data.loc[mol][smi_col])
                                smiles = m.write('can').strip()
                                mw = m.molwt
                                mw = round(mw,3)
                            except:
                                db=db.replace('_',' ')
                                return render_template('insert.html',title=db,err_msg='Invalid SMILES at position %d.'%str(mol))
                            new_entries.append((smiles,'None',mw))
                    else:
                        try:
                            m = pybel.readstring("smiles",data.loc[mol][smi_col])
                            smiles = m.write('can').strip()
                            mw = m.molwt
                            mw = round(mw,3)
                        except:
                            db=db.replace('_',' ')
                            return render_template('insert.html',title=db,err_msg='Invalid SMILES at position %d.'%str(mol))
                        new_entries.append((smiles,data.loc[mol][mol_identifier],mw))

                required_entries = list(set(new_entries) - set(molecules))
                cur.executemany('INSERT INTO MOLECULE(SMILES_str,Molecule_identifier,MW) VALUE(%s,%s,%s)',required_entries)
                print('molecule table populated')
                # populating the credit table
                # todo: figure out how to deal with the credit/publication
                # cur.execute('INSERT INTO %s.Credit(DOI) VALUES(%s)'%db,' ')
                cols=list(conf.properties)
                if smi_col!='':
                    cols.append(smi_col)
                if mol_identifier!='':
                    cols.append(mol_identifier)
                data = data[cols]
                # print(data)
                # populating the values table

                cur.execute('SELECT id,Property_str from Property')
                all_props = cur.fetchall()
                prop_id = dict(map(reversed,all_props))
                cur.execute('SELECT id,Method_name from Model')
                all_models = cur.fetchall()
                model_id = dict(map(reversed,all_models))
                # print(model_id)
                cur.execute('SELECT id,name from Functional')
                all_functionals = cur.fetchall()
                functional_id=dict(map(reversed,all_functionals))
                cur.execute('SELECT id,name from  Basis_set')
                all_basis = cur.fetchall()
                basis_id = dict(map(reversed,all_basis))
                cur.execute('SELECT id,name from Forcefield')
                all_ff = cur.fetchall()
                ff_id = dict(map(reversed,all_ff))
                if smi_col=='':
                    cur.execute("SELECT id,Molecule_identifier from Molecule") 
                    all_mols = cur.fetchall()
                    molecule_id = dict(map(reversed,all_mols))
                    data['molecule_id']=data[mol_identifier].apply(lambda a: molecule_id[a])
                else:
                    cur.execute("SELECT id,SMILES_str from Molecule")
                    all_mols = cur.fetchall()
                    molecule_id = dict(map(reversed,all_mols))
                    data['molecule_id']=data[smi_col].apply(lambda a: molecule_id[pybel.readstring('smi',a).write('can').strip()])
                
                molecule_id = dict(map(reversed,all_mols))
                if mol_identifier!='':
                    data.drop(mol_identifier,1,inplace=True)
                if smi_col!='':
                    data.drop(smi_col,1,inplace=True)
                # print(data.columns)
                # print(data)
                data = data.melt('molecule_id')
                # print(data)
                # print(prop_id)
                data['property_id']=data['variable'].apply(lambda a: prop_id[a])
                # d[df.loc[df['properties'].tolist().index('Density')]['methods']]
                # d[df.loc[df['var'].index.tolist()[1]]['methods']]
                data['model_id']=data['variable'].apply(lambda a: model_id[conf.loc[conf['properties'].tolist().index(a)]['methods']])
                data['functional_id']=data['variable'].apply(lambda a: functional_id[conf.loc[conf['properties'].tolist().index(a)]['functional']])
                data['Basis_id']=data['variable'].apply(lambda a: basis_id[conf.loc[conf['properties'].tolist().index(a)]['basis']])
                data['ff_id']=data['variable'].apply(lambda a: ff_id[conf.loc[conf['properties'].tolist().index(a)]['forcefield']])
                data.drop('variable',1,inplace=True)
                cur.executemany('INSERT INTO VALUE(molecule_id,num_value,property_id,model_id,functional_id,Basis_id,forcefield_id) VALUES(%s,%s,%s,%s,%s,%s,%s)',data.values.tolist())
                print('value table populated')
                con.commit() 
                db=db.replace('_chembddb','')
                db=db.replace('_',' ')
                if host=='':
                    return render_template('insert.html',title=db,success_msg='The database has been successfully populated')
                else:
                    return True

@app.route('/search',methods=['GET','POST'])
def search():
    return render_template('search.html',all_dbs=all_dbs)

@app.route('/search_db<db>',methods=['GET','POST'])
def search_db(db):
    all_dbs=[]
    cur.execute('show databases;')
    all_dbs_tup=cur.fetchall()
    for i in all_dbs_tup:
        if '_chembddb' in i[0]:
            m=i[0]
            all_dbs.append((m[:-9],))
    db=db[1:-1]
    db=db+'_chembddb'
    # print(db)
    cur.execute('USE INFORMATION_SCHEMA')
    result=cur.execute('SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA WHERE SCHEMA_NAME=\'%s\''%db)
    if result == 0:
        error_message='Database does not exist'
        return render_template('Search.html',all_dbs=all_dbs,err_msg=error_message)
    else:
        cur.execute('USE %s'%db)
        db = db.replace('_chembddb','')
        cur.execute('Select * from Property')
        properties=cur.fetchall()
        cur.execute('Select id,method_name from Model')
        results=cur.fetchall()
        cur.execute("Select * from Functional")
        functionals=cur.fetchall()
        cur.execute("Select * from basis_set")
        basis_sets=cur.fetchall()
        cur.execute("Select * from forcefield")
        forcefields=cur.fetchall()
        methods=[]
        global search_results
        is_download=True
        to_order=True
        # print(results)
        # print_l('Preparing the search page...','./')
        for i in results:
            methods.append(i[1])
        # print(methods)
        if request.method == 'POST' and 'search-query' in request.form:
            from_form = request.form
            # print(from_form)
            sql='select value.molecule_id,molecule.SMILES_str,model.method_name,functional.name,basis_set.name,forcefield.name,Property.Property_str, value.num_value from molecule inner join Value on molecule.id=value.molecule_id inner join property on property.id=value.property_id inner join model on model.id=value.model_id inner join functional on functional.id=value.functional_id inner join basis_set on basis_set.id = value.basis_id inner join forcefield on forcefield.id=value.forcefield_id where '
            # sql = 'select value.molecule_id,molecule.SMILES_str,Property.Property_str, value.num_value from molecule inner join Value on molecule.id=value.molecule_id inner join property on property.id=value.property_id where '
            from_form = from_form.to_dict(flat=False)
            prop_counter=0
            keys=[i for i in from_form if '_id' in i]
            min_max_err=False
            min_max_prop=[]
            # print(from_form)
            for k in keys:
                prop_id=int(from_form[k][0])
                from_val=float(from_form[k[:-3]+'_from_val'][0])
                to_val=float(from_form[k[:-3]+'_to_val'][0])
                if from_val > to_val:
                    min_max_err=True
                sql = sql[:sql.rfind('where')+6] + 'molecule_id in (select molecule_id from value where value.property_id={0} and value.num_value>{1} and value.num_value<{2}) and '.format(prop_id,from_val,to_val) + sql[sql.rfind('where')+6:] + 'value.property_id={0} or  '.format(prop_id)
            
            if len(keys)!=0:
                sql=sql[:-5]

            if 'MW' in from_form:
                from_val=float(from_form['MW_from_val'][0])
                to_val=float(from_form['MW_to_val'][0])
                if from_val > to_val:
                    min_max_err=True
                if len(keys)!=0:
                    sql=sql+" and molecule.MW > {} and molecule.MW < {}".format(float(from_form['MW_from_val'][0]),float(from_form['MW_to_val'][0]))
                else:
                    sql=sql+" molecule.MW > {} and molecule.MW < {}".format(float(from_form['MW_from_val'][0]),float(from_form['MW_to_val'][0]))
                    keys.append('MW')

            if 'smiles_search' in from_form:
                if len(keys)==0:
                    sql=sql[:-6]
            
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
            sql=sql+';'
            print(sql)
            cur.execute(sql)
            data1=cur.fetchall() 
            data = pd.DataFrame(list(data1), columns=['Molecule_id','SMILES','Method','Functional','Basis_set','forcefield','Property','Value'])
            data['ID_SMI']=data['Molecule_id'].astype(str)+','+data['SMILES']
            
            data['Property']=data['Property']+'-' +data['Method']+'('+data['Functional']+'/'+data['Basis_set']+')('+data['forcefield']+')'

            data = data[data.columns[-3:]]
            # print(data.head())
            # print(data.columns)
            data=data.pivot_table(index='ID_SMI',columns='Property',values='Value')
            data = data.reset_index()
            # print(data.head())
            # print(data.columns)
            if len(data)>0:
                # data=data['Value'].reset_index()
                data[['ID','SMILES']]=data['ID_SMI'].str.split(',',expand=True)
                columns=['ID','SMILES']
                for i in data.columns[1:-2]:
                    columns.append(i)
                data=data[columns]
                columns=[c.replace('(NA/NA)','') for c in columns]
                columns=[c.replace('(na/na)','') for c in columns]
                columns=[c.replace('(NA)','') for c in columns]
                columns=[c.replace('(na)','') for c in columns]
                try:
                    if 'smiles_search' in from_form:
                        # sql=sql+" and molecule.SMILES_str like \'{}\'".format('%'+from_form['smiles'][0]+'%')
                        smarts = pybel.Smarts(from_form['smiles'][0])
                        # print(smarts)
                        for i in range(len(data)):
                            mol = pybel.readstring("smi",data.loc[i]['SMILES'])
                            smarts.obsmarts.Match(mol.OBMol)
                            if len(smarts.findall(mol))==0:
                                data.drop(i,0,inplace=True)
                    search_results=data
                    search_results.columns=columns
                    if len(data)==0:
                        n_res='Number of results='+ str(len(data))+'\nNo such candidates exist in your database'
                    else:
                        n_res = str(len(data))
                    # print(data)
                    # if 'csv' in from_form:
                    #     search_results.to_csv(from_form['path'][0],index=None)
                except:
                    n_res='Invalid Smarts entered'
                    # print_l('Invalid Smarts entered','.')
                    data=pd.DataFrame()
            else:
                if min_max_err==True:
                    # print_l('Min value entered is > Max value entered for a property','./')
                    n_res = 'Min value entered is > Max value entered for one of the fields above.' 
                    columns=''
                else:               
                    # print_l('No such candidates exist in your database','./')
                    n_res = 'Number of results='+ str(len(data))+'. No such candidates exist in your database'
                    columns=''
            desc=['','']
            for i in data.columns[2:]:
                desc.append('mean={}, std={}, min={}, max={}'.format(data[i].describe()['mean'].round(2),data[i].describe()['std'].round(2),data[i].describe()['min'].round(2),data[i].describe()['max'].round(2)))
            data = tuple(data.itertuples(index=False,name=None))
            is_download=False
            to_order=False
            # print_l('Preparing the results...','./')
            db = db.replace('_chembddb','')
            return render_template('search_db.html', data = data,properties=properties,columns=columns,methods=methods,is_download=is_download,n_res=n_res,basis=basis_sets,functionals=functionals,forcefields=forcefields,all_dbs=all_dbs,title=db,desc=desc)
        elif 'download_csv' in request.form:
            # print(is_download)
            from_form = request.form
            # print(from_form)
            from_form = from_form.to_dict(flat=False)
            # print(search_results)
            desc=['','']
            for i in search_results.columns[2:]:
                desc.append('mean={}, std={}, min={}, max={}'.format(search_results[i].describe()['mean'].round(2),search_results[i].describe()['std'].round(2),search_results[i].describe()['min'].round(2),search_results[i].describe()['max'].round(2)))
            search_results.to_csv('results.csv',index=None)
            msg='Results have been downloaded as results.csv'
            data = tuple(search_results.itertuples(index=False,name=None))
            # print_l('The results have been downloaded in {}.csv ...'.format(from_form['path'][0]),'./')
            return render_template('search_db.html',data = data,properties=properties,columns=search_results.columns,methods=methods,msg=msg,n_res=len(data),functionals=functionals,basis=basis_sets,forcefields=forcefields,is_download=is_download,all_dbs=all_dbs,title=db,desc=desc)
        elif 'orderby_property' in request.form:
            from_form=request.form
            from_form = from_form.to_dict(flat=False)
            # print(search_results)
            if 'ascending' in from_form['select_order']:
                search_results=search_results.sort_values(by=from_form['property_orderby'])
            else:
                search_results=search_results.sort_values(by=from_form['property_orderby'],ascending=False)
            desc=['','']
            data = tuple(search_results.itertuples(index=False,name=None))
            for i in search_results.columns[2:]:
                desc.append('mean={}, std={}, min={}, max={}'.format(search_results[i].describe()['mean'].round(2),search_results[i].describe()['std'].round(2),search_results[i].describe()['min'].round(2),search_results[i].describe()['max'].round(2)))
            n_res=len(data)
            # print_l('The results have been ordered in {} order of {}...'.format(from_form['select_order'],from_form['property_orderby']))
            return render_template('search_db.html',data = data,properties=properties,columns=search_results.columns,methods=methods,n_res=n_res,basis=basis_sets,functionals=functionals,forcefields=forcefields,all_dbs=all_dbs,title=db,desc=desc)        
        else:
            # print(properties)
            return render_template('search_db.html',properties=properties,methods=methods,is_download=is_download,basis=basis_sets, functionals=functionals, forcefields=forcefields,all_dbs=all_dbs,title=db)

@app.route('/molecule-<dbid>',methods=['GET','POST'])
def molecule(dbid):
    
    # cur,conn=connect_mysql()
    import urllib.parse
    
    db=dbid.split('-')[0]
    db=db.replace(' ','_')
    id=dbid.split('-')[1]
    sql = 'SELECT Molecule.id, Molecule.MW, Molecule.SMILES_str,Molecule.Molecule_identifier,Property.Property_str, Property.Unit, Model.method_name,  Functional.name, Basis_set.name,forcefield.name, Value.num_value from Molecule inner join value on Molecule.id=Value.Molecule_id inner join Property on Property.id=VALUE.property_id INNER JOIN Model on Value.model_id=Model.id inner join functional on functional.id=value.functional_id inner join basis_set on basis_set.id=value.basis_id inner join forcefield on forcefield.id=value.forcefield_id where Molecule.id={}'.format(id)
    # print(sql)
    db=db+'_chembddb'
    cur.execute('USE {};'.format(db))
    cur.execute(sql)
    result=cur.fetchall()
    # print(result)
    mol_data=pd.DataFrame(list(result),columns=['ID','MW','SMILES','Identifier','Property','Unit','Method','Functional','Basis_set','Forcefield','Value'])
    mol_data['ALL']=mol_data['ID'].astype(str) +',;'+mol_data['MW'].astype(str)+',;'+mol_data['SMILES']+',;'+mol_data['Identifier']
    mol_data['Property(Unit)']=mol_data['Property']+' ('+mol_data['Unit']+')\n'+'- '+mol_data['Method']+'('+mol_data['Functional']+'/'+mol_data['Basis_set']+')('+mol_data['Forcefield']+')'
    # print(mol_data)
    mol_data=mol_data[['ALL','Property(Unit)','Value']]
    # print(mol_data)
    mol_data=mol_data.pivot(index='ALL',columns='Property(Unit)')
    mol_data=mol_data['Value'].reset_index()
    # print(mol_data)
    mol_data[['ID','MW','SMILES','Identifier']]=mol_data['ALL'].str.split(',;',expand=True)
    url_smi = urllib.parse.quote_plus(mol_data['SMILES'][0])
    # mol_data[['Property(Unit)','Method']]=mol_data['Property(Unit)'].str.split('-',expand=True)
    cols=['ID','MW','SMILES','Identifier']
    for i in mol_data.columns[1:-4]:
        cols.append(i)
    mol_data=mol_data[cols]
    # print(mol_data)
    mol_ob = pybel.readstring("smi",mol_data['SMILES'][0])
    mymol = pybel.readstring("smi", mol_ob.write("can"))
    mymol.make3D(forcefield='mmff94', steps=50)
    # print(os.listdir())
    mymol.write('xyz',app.config['UPLOAD FOLDER']+'/mol_{}.xyz'.format(mol_data['ID'][0]),overwrite=True)
    # mymol.write('xyz',app.config['UPLOAD FOLDER']+'mol.xyz',overwrite=True)
    # f = open('./mol_{}.xyz'.format(mol_data['ID'][0]))
    cols=[c.replace('(na/na)','') for c in cols]
    cols=[c.replace('(na)','') for c in cols]
    cols=[c.replace('(na)','') for c in cols]
    mol_data = tuple(mol_data.itertuples(index=False,name=None))
    mol_data = (tuple(cols[:]),)+mol_data
    mol_data=tuple(zip(*mol_data))
    db=db.replace('_',' ')
    return render_template('molecule.html',mol_data=mol_data,columns=cols,title=db,all_dbs=all_dbs,url_smi=url_smi)

@app.route('/delete',methods=['GET','POST'])
def delete(host='',user='',pw='',db=''):
    if db !='':
        a,b = connect_mysql(host=host,user=user,pw=pw)
        cur.execute('drop database %s;'%db)
        return True
    else:
        cur.execute('show databases;')
        all_dbs=[]
        all_dbs_tup=cur.fetchall()
        for i in all_dbs_tup:
            if '_chembddb' in i[0]:
                m=i[0]
                all_dbs.append((m[:-9],))
        details=request.form
        details=details.to_dict(flat=True)
        # print(details)
        if 'dbname' in details:
            # print(basis_sets)
            dbname=details['dbname']+'_chembddb'
            cur.execute('use {};'.format(dbname))
            cur.execute('Select * from Property')
            properties=cur.fetchall()
            cur.execute('Select id,method_name from Model')
            results=cur.fetchall()
            cur.execute("Select * from Functional")
            functionals=cur.fetchall()
            cur.execute("Select * from basis_set")
            basis_sets=cur.fetchall()
            cur.execute("Select * from forcefield")
            forcefields=cur.fetchall()
            methods=[]
            for i in results:
                methods.append(i[1])
        if 'submit' in details:
            return render_template('delete.html',data=True,dbname=dbname,properties=properties,methods=methods,functionals=functionals,basis=basis_sets,forcefields=forcefields,all_dbs=all_dbs)
        elif 'search-query' in details:

            if 'exampleRadios' not in details:
                cur.execute('drop database {}'.format(details['dbname']))
                # print('done')
                cur.execute('show databases;')
                all_dbs_tup=cur.fetchall()
                all_dbs=[]
                for i in all_dbs_tup:
                    if '_chembddb' in i[0]:
                        m=i[0]
                        all_dbs.append((m[:-9],))
                return render_template('delete.html',data=True,properties=properties,methods=methods,functionals=functionals,basis=basis_sets,forcefields=forcefields,all_dbs=all_dbs,success_msg='database {} deleted'.format(details['dbname']))

            elif details['exampleRadios']=='option1':
                if details['MW']!='':
                    mw_from=details['MW_from_val']
                    mw_to=details['MW_to_val']

            else:
                cur.execute('USE {};'.format(details['dbname']))
                keys=[i for i in details if '_id' in i]
                # find molecule id if smiles_search in details
                # find MW from and to
                
                for k in keys:
                    prop_id=int(details[k][0])
                    # print(prop_id)
                    if details[k[:-3]+'_from_val']=='' and details[k[:-3]+'_to_val']=='':
                        # print('here')
                        sql='DELETE FROM Property WHERE id={};'.format(prop_id)
                        # print(sql)
                        cur.execute(sql)
                        cur.execute('use {};'.format(dbname))
                        cur.execute('Select * from Property')
                        properties=cur.fetchall()
                        cur.execute('Select id,method_name from Model')
                        results=cur.fetchall()
                        cur.execute("Select * from Functional")
                        functionals=cur.fetchall()
                        cur.execute("Select * from basis_set")
                        basis_sets=cur.fetchall()
                        cur.execute("Select * from forcefield")
                        forcefields=cur.fetchall()
                        methods=[]
                        for i in results:
                            methods.append(i[1])
                    else:
                        from_val=float(details[k[:-3]+'_from_val'])
                        # print(from_val)
                        to_val=float(details[k[:-3]+'_to_val'])
                        # print(to_val)
                        if from_val > to_val:
                            return render_template('delete.html',data=True, dbname=details['dbname'].replace('_chembddb',''),properties=properties,methods=methods,functionals=functionals,basis=basis_sets,forcefields=forcefields,all_dbs=all_dbs,err_msg='Minimum value for one of the properties is greater than the maximum value for it.')
                        else:
                            sql='DELETE FROM Value WHERE property_id={} and num_value > {} and num_value < {};'.format(prop_id,from_val,to_val)
                            # print(sql)
                            cur.execute(sql)
            return render_template('delete.html',data=True, dbname=details['dbname'].replace('_chembddb',''),properties=properties,methods=methods,functionals=functionals,basis=basis_sets,forcefields=forcefields,all_dbs=all_dbs,success_msg='Deleted from database {}.'.format(details['dbname'].replace('_chembddb','')))
        else:
            return render_template('delete.html',all_dbs=all_dbs)

@app.route('/uploads/<filename>')
def uploaded_file(filename):
    return send_from_directory(app.config['UPLOAD FOLDER'],filename)

def run_config():
    print('Open http://localhost:5000/connect')
    app.run(debug=True)