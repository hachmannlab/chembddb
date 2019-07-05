from flask import Flask, render_template, url_for, request
import pymysql
import os
import sys
import pandas as pd
from copy import deepcopy
import pybel

app = Flask(__name__)
def connect_mysql(host,user,pw):
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
    con = pymysql.connect(host = host, user=user, password = pw)
    cur = con.cursor()
    return cur, con

def setup(cur,db):
    """Function to create the database with all the tables (Molecule, Credit/Publication, Property, Model, Topology, Value), with all their rows and their acceptable value types. 

    Parameters
    ----------
    cursor: pymysql cursor object
        pointer to the MySQL database; used to execute all SQL queries
    db: str
        name of the database that needs to be set up
    Returns
    -------

    """

    # con = pymysql.connect(host = host, user = user, password = pw)
    # cur = con.cursor()
    # cur.execute('DROP DATABASE IF EXISTS %s;'%db)
    cur.execute('USE INFORMATION_SCHEMA')
    result = cur.execute('SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA WHERE SCHEMA_NAME=\'%s\''%db)
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

        cur.execute('ALTER TABLE `%s`.`Value` ADD CONSTRAINT `Value_fk0` FOREIGN KEY (`model_id`) REFERENCES `Model`(`id`);'%db)
        cur.execute('ALTER TABLE `%s`.`Value` ADD CONSTRAINT `Value_fk1` FOREIGN KEY (`property_id`) REFERENCES `Property`(`id`);'%db)
        cur.execute('ALTER TABLE `%s`.`Value` ADD CONSTRAINT `Value_fk2` FOREIGN KEY (`molecule_id`) REFERENCES `Molecule`(`id`);'%db)
        cur.execute('ALTER TABLE `%s`.`Value` ADD CONSTRAINT `Value_fk3` FOREIGN KEY (`functional_id`) REFERENCES `Functional`(`id`);'%db)
        cur.execute('ALTER TABLE `%s`.`Value` ADD CONSTRAINT `Value_fk4` FOREIGN KEY (`basis_id`) REFERENCES `Basis_set`(`id`);'%db)
        cur.execute('ALTER TABLE `%s`.`Value` ADD CONSTRAINT `Value_fk5` FOREIGN KEY (`forcefield_id`) REFERENCES `Forcefield`(`id`);'%db)
        # cur.ex ecute('ALTER TABLE `%s`.`Value` ADD CONSTRAINT `Value_fk3` FOREIGN KEY (`credit_id`) REFERENCES `Credit`(`id`);'%db)
        print_l('--The database %s has been setup--'%db,'.')
    else:
        tmp_str = 'ERROR: Database with name %s already exists'%db
        print_le(tmp_str,'.','Change the database name in the config file. Aborting due to improper naming of the database.')

def insert(con, cur,db, inputs):
    """Set docstring here.

    Parameters
    ----------
    con: pymysql connection object
        helps connect to the MySQL server
    cur:  pymysql cursor object
        pointer to the MySQL database; used to execute all SQL queries
    db: str
        database name
    inputs: dict
        contains all the task options input by the user

    Returns
    -------

    """
    # if the user only runs the insert task (without setup), the following lines check if the database exists in the local machine
    cur.execute('USE INFORMATION_SCHEMA')
    result = cur.execute('SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA WHERE SCHEMA_NAME=\'%s\''%db)
    if result == 0:
        tmp_str = 'ERROR: Database %s does not exist.'%db
        print_le(tmp_str,'.','Setup the database or change the name. Aborting due to improper naming of the database.')
    else:
        cur.execute('USE %s;'%db)
        path = inputs['Path to csv']
        data = pd.read_csv(path)
        # n_u=len(inputs['Units list'])
        # n_p=len(inputs['Properties list'])
        # n_m=len(inputs['Methods list'])
        # n_f=
        if len(inputs['Units list'])!=len(inputs['Properties list'])!=len(inputs['Methods list'])!=len(inputs['Functional list'])!=len(inputs['Basis set list'])!=len(inputs['Forcefield list']):
            tmp_str = 'ERROR: Each property does not have a unit/method specified for it.'
            print_le(tmp_str,'.','List a unit/method for each property. Aborting due to different number of properties and units/methods being provided.')

        # loop throught he CSV file, check if the smiles value is in the table, if yes, fetch the corresponding id, same goes for property, same goes for method for that property, if it does not exist, fetch the last id and create a new entry

        # populating and property table
        
        entered_list=[]
        cur.execute("SELECT Property_str from Property")
        properties = cur.fetchall()
        # print(properties)
        for prop, units in zip(inputs['Properties list'],inputs['Units list']):
            # print(prop)
            # print(properties)
            if any(prop in i for i in properties) or prop in entered_list:
                pass
            else:
                entered_list.append(prop)
                print(entered_list)
                cur.execute("INSERT INTO Property(Property_str,Unit) VALUES(%s,%s)",[prop,units])

        print_l('--Property table has been populated--','./')
        # populating the model table
        cur.execute("SELECT Method_name from Model")
        models = cur.fetchall()
        entered_list=[]
        for method in inputs['Methods list']:
            if any(method in i for i in models) or method in entered_list:
                pass
            else:
                entered_list.append(method)
                cur.execute("INSERT INTO Model(Method_name) VALUES(%s)",[method])
        print_l('--Model table has been populated--','./')

        # populating the functional table
        cur.execute("SELECT name FROM functional")
        functionals = cur.fetchall()
        entered_list = []
        for func in inputs['Functional list']:
            if any(func in i for i in functionals) or func in entered_list:
                pass
            else:
                entered_list.append(func)
                cur.execute("INSERT INTO Functional(name) VALUES(%s)",[func])
        print_l('--Functional table has been populated--','./')

        # populating the basis_set table
        cur.execute("SELECT name FROM Basis_set")
        basis_sets = cur.fetchall()
        entered_list = []
        for basis in inputs['Basis set list']:
            if any(basis in i for i in basis_sets) or basis in entered_list:
                pass
            else:
                entered_list.append(basis)
                cur.execute("INSERT INTO Basis_set(name) VALUES(%s)",[basis])
        print_l('--Basis_sets table has been populated--','./')

        # populating the forcefield table
        cur.execute("SELECT name FROM Forcefield")
        forcefields = cur.fetchall()
        entered_list = []
        for ff in inputs['Force Field list']:
            if any(ff in i for i in forcefields) or ff in entered_list or ff=='None':
                pass
            else:
                entered_list.append(basis)
                cur.execute("INSERT INTO Forcefield(name) VALUES(%s)",[ff])
        print_l('--Forcefield table has been populated--','./')

        # populating the molecule table
        cur.execute("SELECT SMILES_str,Molecule_identifier,MW from Molecule")
        molecules = cur.fetchall()
        new_entries=[]
        # print(data[data.filter(regex='smiles|inchi').columns])
        # check last id for molecule, add molecule index, melt dataframe, add property and method index using lambda and dictionary
        # check for molecule identifier (name, IUPAC, )
        
        try:
            smi_ind=list(map(str.lower,inputs['Molecule identifier list'])).index('smiles')
        except:
            smi_ind=-1
        n_identifiers=len(inputs['Molecule identifier list'])
        for mol in range(len(data)):
            # data[data.filter(regex='smiles').columns]
            
            if smi_ind == -1:
                if data.loc[mol][inputs['Molecule identifier list'][0]] not in new_entries:
                    new_entries.append(('None',data.loc[mol][inputs['Molecule identifier list'][0]],'None'))
                else:
                    pass
            else:
                try:
                    m = pybel.readstring("smiles",data.loc[mol][inputs['Molecule identifier list'][smi_ind]])
                    smiles = m.write('can').strip()
                    mw = m.molwt
                    mw = round(mw,3)
                except:
                    tmp_str = 'Error: The SMILES provided in data file is invalid.'
                    print_le(tmp_str, './',"Please provide correct SMILES. Aborting due to wrong molecule description.")
            
                if n_identifiers>1:
                    if smiles in new_entries or data.loc[mol][inputs['Molecule identifier list'][n_identifiers-smi_ind-1]] not in new_entries:
                        pass
                    else:
                        new_entries.append((smiles,data.loc[mol][inputs['Molecule identifier list'][n_identifiers-smi_ind-1]],'None'))
                else:
                    if smiles in new_entries:
                        pass
                    else:
                        new_entries.append((smiles,'None',mw))
        # print(new_entries[0])
        # print(molecules[0])
        required_entries = list(set(new_entries) - set(molecules))
        # print(len(required_entries))
        # print_le('abc','./','abc')
        cur.executemany('INSERT INTO MOLECULE(SMILES_str,Molecule_identifier,MW) VALUE(%s,%s,%s)',required_entries)
        print_l('--Molecule table has been populated--','./')

        # populating the credit table
        # todo: figure out how to deal with the credit/publication
        # cur.execute('INSERT INTO %s.Credit(DOI) VALUES(%s)'%db,' ')


        # populating the values table
        # fetching all entries
        cols = deepcopy(inputs['Properties list'])
        cols = cols + inputs['Molecule identifier list']
        # print(data)
        # print(cols)
        data = data[cols]

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
        if smi_ind != -1:
            cur.execute("SELECT id,SMILES_str from Molecule")
        else:
            cur.execute("SELECT id,Molecule_identifier from Molecule")

        all_mols = cur.fetchall()
        molecule_id = dict(map(reversed,all_mols))
        # print(functional_id)
        # print(inputs['Functional list'][inputs['Methods list'].index()])
        # print(properties)
        if smi_ind !=-1:
            data['molecule_id']=data[data.columns[-1]].apply(lambda a: molecule_id[pybel.readstring('smi',a).write('can').strip()])
        else:
            data['molecule_id']=data[data.columns[-1]].apply(lambda a: molecule_id[a])

        data.drop(data[inputs['Molecule identifier list']],1,inplace=True)
        data = data.melt('molecule_id')
        # print(data)
        data['property_id']=data['variable'].apply(lambda a: prop_id[a])
        # print(inputs['Properties list'])
        data['model_id']=data['variable'].apply(lambda a: model_id[inputs['Methods list'][inputs['Properties list'].index(a)]])
        data['functional_id']=data['variable'].apply(lambda a: functional_id[inputs['Functional list'][inputs['Properties list'].index(a)]])
        data['Basis_id']=data['variable'].apply(lambda a: basis_id[inputs['Basis set list'][inputs['Properties list'].index(a)]])
        # print(ff_id)
        # print(inputs['Force Field list'])
        data['ff_id']=data['variable'].apply(lambda a: ff_id[inputs['Force Field list'][inputs['Properties list'].index(a)]])
        data.drop('variable',1,inplace=True)
        cur.executemany('INSERT INTO VALUE(molecule_id,num_value,property_id,model_id,functional_id,Basis_id,forcefield_id) VALUES(%s,%s,%s,%s,%s,%s,%s)',data.values.tolist())
        con.commit()
        cur.close()

def get_options(config):
    """ 
    Function to read options provided in the config file.
        
    Parameters
    ----------
    config_file: file handle
    Returns
    -------
    options_dict: dict
        dictionary of generation rules provided by the user. if the user provides default values for any rule, it is not added to the dictionary.
    task_options: list
        list of other input arguments to the library generator
    """
    options_list, task_options = [], {'insert':{},'search':{}}
    task_flag = False
    task=''
    output_dir='.'
    for i,line in enumerate(config):
        if i == 0:
            print_l('--CONFIG FILE--', output_dir)
            print_l('-Input_options-', output_dir)
            continue
        elif '==' in line:
            words = line.split('==')
            value = words[1].strip()
        
            if value == 'None':
                tmp_str = "ERROR: Option not provided for "+words[0].strip()
                print_le(tmp_str,output_dir,"Aborting due to incomplete optipns.")

            else:
                if ',' in value:
                    value = ''.join(value.split())
                    options_list.append(value.split(','))
                else:
                    options_list.append(value)
                print_l(line[:-1], output_dir)
        
        elif '##' in line:
            if line.split()[1] in options_list[4]:
                count = 0
                task = line.split()[1]
                print_l(line[:-1], output_dir)
                task_flag = True
            else:
                task_flag = False
                continue

        if task_flag == True:
            if 'end' in line:
                # print(count)
                if count <=0:
                    tmp_str= "ERROR: No task options provided"
                    print_le(tmp_str,output_dir,'Aborting due to incomplete task options.')
                else:
                    print_l(line[:-1], output_dir)
            elif '::' in line:
                words = line.split('::')
                value = words[1].strip()
                # print(value)
                if 'None' in value:
                    if options_list[4]=='insert':
                        # if 'csv' not in task_options['insert']['Input mode'] and 'Path' in words[0] and value == 'None':
                            # pass
                        # else:
                        tmp_str = 'ERROR: The option for the subtask \"'+words[0].strip()+'\" cannot be None.'
                        print_le(tmp_str,'.','Provide options for this subtask. Aborting due to option not provided.')
                    count = count - 1
                elif 'False' in value:
                    count = count - 1
                else:
                    count = count + 1
                if 'list' in words[0]:
                    value = ''.join(value.split())
                    task_options[task][words[0].split('.')[1].strip()]=value.split(',')
                else:
                    task_options[task][words[0].split('.')[1].strip()]=value
                print_l(line[:-1], output_dir)   
        else:
            continue

    return options_list,task_options

def print_l(sentence, output_dir):
    """Print to logfile.
    Parameters
    ----------
    sentence: str
        string to be printed to logfile
    Returns
    -------
    """ 
    # os.makedirs(output_dir)

    logfile = open(os.path.join(output_dir+'/logfile.txt'),'a')
    print(sentence)
    logfile.write(str(sentence)+"\n")

def print_le(sentence, output_dir, msg="Aborting the run"):
    """Print to both error file and logfile and then exit code.
    Parameters
    ----------
    sentence: str
        string to be printed to error file and logfile
    Returns
    -------
    """

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    logfile = open(os.path.join(output_dir+'/logfile.txt'),'a')
    error_file = open(os.path.join(output_dir+'/error_file.txt'),'a')
    print(sentence)
    logfile.write(sentence+"\n")
    error_file.write(sentence+"\n")
    sys.exit(msg)

@app.route('/search',methods=['GET','POST'])
def search():
    cursor.execute('USE %s'%options[3])
    cursor.execute('Select * from Property')
    properties=cursor.fetchall()
    cursor.execute('Select id,method_name from Model')
    results=cursor.fetchall()
    cursor.execute("Select * from Functional")
    functionals=cursor.fetchall()
    cursor.execute("Select * from basis_set")
    basis_sets=cursor.fetchall()
    cursor.execute("Select * from forcefield")
    forcefields=cursor.fetchall()
    methods=[]
    global search_results
    is_download=True
    to_order=True
    # print(results)
    print_l('Preparing the search page...','./')
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

        for k in keys:
            prop_id=int(from_form[k][0])
            from_val=float(from_form[k[:-3]+'_from_val'][0])
            to_val=float(from_form[k[:-3]+'_to_val'][0])
            sql = sql[:sql.rfind('where')+6] + 'molecule_id in (select molecule_id from value where value.property_id={0} and value.num_value>{1} and value.num_value<{2}) and '.format(prop_id,from_val,to_val) + sql[sql.rfind('where')+6:]
        
        if len(keys)!=0:
            sql=sql[:-5]

        if 'MW' in from_form:
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
        cursor.execute(sql)
        data1=cursor.fetchall() 
        data = pd.DataFrame(list(data1), columns=['Molecule_id','SMILES','Method','Functional','Basis_set','forcefield','Property','Value'])
        # print(data)
        # property_values={}
        # for prop in properties:
        #     # print(prop)
        #     p=[]
        #     for i in range(len(data)):
        #         if data.loc[i]['Property']==prop[1]:
        #             p.append(data.loc[i]['Value'])
        #     property_values[prop[1]]=p

        # smi=data['SMILES'][:int(len(data)/len(properties))]
        # data=pd.DataFrame(smi,columns=['SMILES'])
        # for k,v in property_values.items():
        #     data[k]=v
        data['ID_SMI']=data['Molecule_id'].astype(str)+','+data['SMILES']
        data['Property']=data['Property']+'-' +data['Method']+'('+data['Functional']+'/'+data['Basis_set']+')('+data['forcefield']+')'
        print(data)
        data=data.pivot(index='ID_SMI',columns='Property')
        data=data['Value'].reset_index()
        data[['ID','SMILES']]=data['ID_SMI'].str.split(',',expand=True)
        columns=['ID','SMILES']
        for i in data.columns[1:-2]:
            columns.append(i)
        data=data[columns]
        columns=[c.replace('(NA/NA)','') for c in columns]
        columns=[c.replace('(NA)','') for c in columns]
        if 'smiles_search' in from_form:
            # sql=sql+" and molecule.SMILES_str like \'{}\'".format('%'+from_form['smiles'][0]+'%')
            smarts = pybel.Smarts(from_form['smiles'][0])
            for i in range(len(data)):
                mol = pybel.readstring("smi",data.loc[i]['SMILES'])
                smarts.obsmarts.Match(mol.OBMol)
                if len(smarts.findall(mol))==0:
                    data.drop(i,0,inplace=True)
        search_results=data
        search_results.columns=columns
        # print(data)
        # if 'csv' in from_form:
        #     search_results.to_csv(from_form['path'][0],index=None)

        data = tuple(data.itertuples(index=False,name=None))
        is_download=False
        to_order=False
        n_res = len(data)
        print_l('Preparing the results...','./')
        return render_template('index.html', data = data,properties=properties,columns=columns,title=options[3],methods=methods,is_download=is_download,n_res=n_res,basis=basis_sets,functionals=functionals,forcefields=forcefields)
    elif 'download_csv' in request.form:
        # print(is_download)
        from_form = request.form
        # print(from_form)
        from_form = from_form.to_dict(flat=False)
        # print(search_results)
        search_results.to_csv(from_form['path'][0]+'.csv',index=None)
        msg='Results have been downloaded as {}.csv'.format(from_form['path'][0])
        data = tuple(search_results.itertuples(index=False,name=None))
        print_l('The results have been downloaded in {}.csv ...'.format(from_form['path'][0]),'./')
        return render_template('index.html',data = data,properties=properties,columns=search_results.columns,title=options[3],methods=methods,msg=msg,is_download=is_download)
    elif 'orderby_property' in request.form:
        from_form=request.form
        from_form = from_form.to_dict(flat=False)
        print(search_results)
        if 'ascending' in from_form['select_order']:
            search_results=search_results.sort_values(by=from_form['property_orderby'])
        else:
            search_results=search_results.sort_values(by=from_form['property_orderby'],ascending=False)
        data = tuple(search_results.itertuples(index=False,name=None))
        n_res=len(data)
        # print_l('The results have been ordered in {} order of {}...'.format(from_form['select_order'],from_form['property_orderby']))
        return render_template('index.html',data = data,properties=properties,columns=search_results.columns,title=options[3],methods=methods,n_res=n_res,basis=basis_sets,functionals=functionals,forcefields=forcefields)        
    else:
        # print(properties)
        return render_template('index.html',properties=properties,methods=methods,title=options[3],is_download=is_download,basis=basis_sets, functionals=functionals, forcefields=forcefields)

@app.route('/molecule-<id>',methods=['GET','POST'])
def molecule(id):
    sql = 'SELECT Molecule.id, Molecule.MW, Molecule.SMILES_str,Molecule.Molecule_identifier,Property.Property_str, Property.Unit, Model.method_name,  Functional.name, Basis_set.name,forcefield.name, Value.num_value from Molecule inner join value on Molecule.id=Value.Molecule_id inner join Property on Property.id=VALUE.property_id INNER JOIN Model on Value.model_id=Model.id inner join functional on functional.id=value.functional_id inner join basis_set on basis_set.id=value.basis_id inner join forcefield on forcefield.id=value.forcefield_id where Molecule.id={}'.format(id)
    # print(sql)
    cursor.execute('USE {};'.format(options[3]))
    cursor.execute(sql)
    result=cursor.fetchall()
    # print(result)
    mol_data=pd.DataFrame(list(result),columns=['ID','MW','SMILES','Identifier','Property','Unit','Method','Functional','Basis_set','Forcefield','Value'])
    mol_data['ALL']=mol_data['ID'].astype(str) +',;'+mol_data['MW'].astype(str)+',;'+mol_data['SMILES']+',;'+mol_data['Identifier']
    mol_data['Property(Unit)']=mol_data['Property']+' ('+mol_data['Unit']+')\n'+'- '+mol_data['Method']+'('+mol_data['Functional']+'/'+mol_data['Basis_set']+')('+mol_data['Forcefield']+')'
    print(mol_data)
    mol_data=mol_data[['ALL','Property(Unit)','Value']]
    print(mol_data)
    mol_data=mol_data.pivot(index='ALL',columns='Property(Unit)')
    mol_data=mol_data['Value'].reset_index()
    print(mol_data)
    mol_data[['ID','MW','SMILES','Identifier']]=mol_data['ALL'].str.split(',;',expand=True)
    
    # mol_data[['Property(Unit)','Method']]=mol_data['Property(Unit)'].str.split('-',expand=True)
    cols=['ID','MW','SMILES','Identifier']
    for i in mol_data.columns[1:-4]:
        cols.append(i)
    mol_data=mol_data[cols]
    # print(mol_data)
    mol_ob = pybel.readstring("smi",mol_data['SMILES'][0])
    mymol = pybel.readstring("smi", mol_ob.write("can"))
    mymol.make3D(forcefield='mmff94', steps=50)
    mymol.write('xyz', './static/xyz/mol_{}.xyz'.format(mol_data['ID'][0]),overwrite=True)
    cols=[c.replace('(NA/NA)','') for c in cols]
    cols=[c.replace('(NA)','') for c in cols]
    mol_data = tuple(mol_data.itertuples(index=False,name=None))
    mol_data = (tuple(cols[:]),)+mol_data
    mol_data=tuple(zip(*mol_data))

    return render_template('molecule.html',mol_data=mol_data,columns=cols,title=options[3])
def run_config(input_file):
    f = open(input_file)
    global cursor, connection, options
    options,subtasks=get_options(f)
    f.close()
    cursor,connection = connect_mysql(host=options[0],user=options[1],pw=options[2])
    tasks = options[4]
    if 'setup' in tasks:
        setup(cursor,options[3])
        if 'search' in tasks and len(tasks) == 2:
            tmp_str = 'Error: Database is empty.'
            print_le(tmp_str,'.','Insert values into the database and then search. Aborting due to empty database')
    
    if 'insert' in tasks:
        insert(connection, cursor, options[3],subtasks['insert'])
    
    if 'search' in tasks:
        app.run(debug=True)


