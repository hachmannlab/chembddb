from flask import Flask, render_template, url_for, request
import pymysql
import os
import sys
import pandas as pd
from copy import deepcopy
import pybel
from flask import send_from_directory
import numpy as np

app = Flask(__name__)
# upload_directory='/Users/adi/chembddb-1/'
upload_directory=os.getcwd()
app.config['UPLOAD FOLDER']=upload_directory


def connect_mysql(cred={'host':'','username':'','password':''}):
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
    if os.path.isfile(app.config['UPLOAD FOLDER']+'/credentials.dat'):
        with open(app.config['UPLOAD FOLDER']+'/credentials.dat') as file:
            lines=file.readlines()
        cred_from_file={}
        for line in lines:
            cred_from_file[line.split(':')[0].strip()]=line.split(':')[1].strip()
    # print(cred_from_file)
    
    for key in cred.keys():
        if cred[key]=='':
            cred[key]=cred_from_file[key]
    # print(cred)
    con = pymysql.connect(host = cred['host'], user=cred['username'], password = cred['password'])
    cur = con.cursor()
    return cur, con

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

@app.route('/setup',methods=['GET','POST'])
def setup():
    if request.method=='POST' and 'setup' in request.form:
        cred=request.form
        cred = cred.to_dict(flat=False)
        # print(cred)
        # storing the credentials in a credentials.dat
        if 'first-time' in cred:
            with open(app.config['UPLOAD FOLDER']+'/credentials.dat','w') as file:
                file.write('host:'+cred['host'][0]+'\n')
                file.write('username:'+cred['username'][0]+'\n')
                file.write('password:'+cred['password'][0])
            cur,conn = connect_mysql({'host':cred['host'][0],'username':cred['username'][0],'password':cred['password'][0]})
        else:
            cur,conn = connect_mysql()
        # connecting to mysql
        

        db=cred['dbname'][0]
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

    return render_template('setup.html')

@app.route('/insert<db>',methods=['GET','POST'])
def insert(db):
    db=db[1:-1]
    mi_cols=[]
    cur,conn=connect_mysql()
    cur.execute('USE INFORMATION_SCHEMA')
    result = cur.execute('SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA WHERE SCHEMA_NAME=\'%s\''%db)
    if result == 0:
        error_message='Database does not exist'
        return render_template('setup.html',title=db,err_msg=error_message)
    else:
        if request.method=='POST':
            config_options=request.form
            print(config_options)
            config_options=config_options.to_dict(flat=False)
            files=request.files
            files=files.to_dict(flat=False)
            print(files)
            conf_file=files['config_file'][0]
            data_file=files['data_file'][0]
            smi_col=config_options['smiles'][0]
            mol_identifier=config_options['molecule_identifier'][0]
            print(mol_identifier)

            if conf_file.filename=='' or conf_file.filename.rsplit('.',1)[1]!='csv':
                return render_template('insert.html',title=db,err_msg='No config file provided or incorrect file format. (csv required)')
            elif data_file.filename=='' or data_file.filename.rsplit('.',1)[1]!='csv':
                return render_template('insert.html',title=db,err_msg='No data file provided or incorrect file format. (csv required)')
            elif smi_col=='' and mol_identifier=='':
                return render_template('insert.html',title=db,err_msg='No molecule identifiers provided.')
            else:
                conf=pd.read_csv(conf_file)
                print(conf)
                data=pd.read_csv(data_file)
                conf.replace(np.nan,'na',inplace=True)
                all_prop=True
                all_mols=True
                for prop in conf['properties']:
                    if prop not in data.columns:
                        all_prop=False
                if all_prop==False:
                    return render_template('insert.html',title=db,err_msg='Property(s) in config file do not exist in data.')
                else:
                    if smi_col!='' and smi_col not in data.columns:
                        all_mols=False
                    if mol_identifier!='' and mol_identifier not in data.columns:
                        all_mols=False
                    if all_mols==False:
                        return render_template('insert.html',title=db,err_msg='Identifier(s) listed do not exist in data.')
                    else:
                        cur.execute('USE %s;'%db)

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
                                        return render_template('insert.html',title=db,err_msg='Invalid SMILES at position %d.'%str(mol))
                                    new_entries.append((smiles,'None',mw))
                            else:
                                try:
                                    m = pybel.readstring("smiles",data.loc[mol][smi_col])
                                    smiles = m.write('can').strip()
                                    mw = m.molwt
                                    mw = round(mw,3)
                                except:
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
                        # populating the values table

                        cur.execute('SELECT id,Property_str from Property')
                        all_props = cur.fetchall()
                        prop_id = dict(map(reversed,all_props))
                        cur.execute('SELECT id,Method_name from Model')
                        all_models = cur.fetchall()
                        model_id = dict(map(reversed,all_models))
                        print(model_id)
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
                        data.drop(cols[-2:],1,inplace=True)
                        data = data.melt('molecule_id')
                        print(data)
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
                        conn.commit() 
                        
                        return render_template('insert.html',title=db,success_msg='The database has been successfully populated')
        return render_template('insert.html',title=db)

@app.route('/search<db>',methods=['GET','POST'])
def search(db):
    db=db[1:-1]
    cursor,con=connect_mysql()
    cursor.execute('USE INFORMATION_SCHEMA')
    result = cursor.execute('SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA WHERE SCHEMA_NAME=\'%s\''%db)
    if result == 0:
        error_message='Database does not exist'
        return render_template('setup.html',title=db,err_msg=error_message)
    else:
        cursor.execute('USE %s'%db)
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
            min_max_err=False
            min_max_prop=[]
            for k in keys:
                prop_id=int(from_form[k][0])
                from_val=float(from_form[k[:-3]+'_from_val'][0])
                to_val=float(from_form[k[:-3]+'_to_val'][0])
                if from_val > to_val:
                    min_max_err=True
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
            # print(data)
            data=data.pivot(index='ID_SMI',columns='Property')
            if len(data)>0:
                data=data['Value'].reset_index()
                data[['ID','SMILES']]=data['ID_SMI'].str.split(',',expand=True)
                columns=['ID','SMILES']
                for i in data.columns[1:-2]:
                    columns.append(i)
                data=data[columns]
                columns=[c.replace('(NA/NA)','') for c in columns]
                columns=[c.replace('(NA)','') for c in columns]
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
                    print_l('Invalid Smarts entered','.')
                    data=pd.DataFrame()
            else:
                if min_max_err==True:
                    print_l('Min value entered is > Max value entered for a property','./')
                    n_res = 'Min value entered is > Max value entered for a property' 
                    columns=''
                else:               
                    print_l('No such candidates exist in your database','./')
                    n_res = 'Number of results='+ str(len(data))+'\nNo such candidates exist in your database'
                    columns=''
            data = tuple(data.itertuples(index=False,name=None))
            is_download=False
            to_order=False
            
            print_l('Preparing the results...','./')
            return render_template('index.html', data = data,properties=properties,columns=columns,title=db,methods=methods,is_download=is_download,n_res=n_res,basis=basis_sets,functionals=functionals,forcefields=forcefields)
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
            return render_template('index.html',data = data,properties=properties,columns=search_results.columns,title=db,methods=methods,msg=msg,is_download=is_download)
        elif 'orderby_property' in request.form:
            from_form=request.form
            from_form = from_form.to_dict(flat=False)
            # print(search_results)
            if 'ascending' in from_form['select_order']:
                search_results=search_results.sort_values(by=from_form['property_orderby'])
            else:
                search_results=search_results.sort_values(by=from_form['property_orderby'],ascending=False)
            data = tuple(search_results.itertuples(index=False,name=None))
            n_res=len(data)
            # print_l('The results have been ordered in {} order of {}...'.format(from_form['select_order'],from_form['property_orderby']))
            return render_template('index.html',data = data,properties=properties,columns=search_results.columns,title=db,methods=methods,n_res=n_res,basis=basis_sets,functionals=functionals,forcefields=forcefields)        
        else:
            # print(properties)
            return render_template('index.html',properties=properties,methods=methods,title=db,is_download=is_download,basis=basis_sets, functionals=functionals, forcefields=forcefields)

@app.route('/molecule-<dbid>',methods=['GET','POST'])
def molecule(dbid):
    cursor,conn=connect_mysql()
    
    db=dbid.split('-')[0]
    id=dbid.split('-')[1]
    sql = 'SELECT Molecule.id, Molecule.MW, Molecule.SMILES_str,Molecule.Molecule_identifier,Property.Property_str, Property.Unit, Model.method_name,  Functional.name, Basis_set.name,forcefield.name, Value.num_value from Molecule inner join value on Molecule.id=Value.Molecule_id inner join Property on Property.id=VALUE.property_id INNER JOIN Model on Value.model_id=Model.id inner join functional on functional.id=value.functional_id inner join basis_set on basis_set.id=value.basis_id inner join forcefield on forcefield.id=value.forcefield_id where Molecule.id={}'.format(id)
    # print(sql)
    cursor.execute('USE {};'.format(db))
    cursor.execute(sql)
    result=cursor.fetchall()
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
    cols=[c.replace('(NA/NA)','') for c in cols]
    cols=[c.replace('(NA)','') for c in cols]
    mol_data = tuple(mol_data.itertuples(index=False,name=None))
    mol_data = (tuple(cols[:]),)+mol_data
    mol_data=tuple(zip(*mol_data))
    return render_template('molecule.html',mol_data=mol_data,columns=cols,title=db)

@app.route('/uploads/<filename>')
def uploaded_file(filename):
    return send_from_directory(app.config['UPLOAD FOLDER'],filename)

def run_config(option):
    app.run(debug=True)
    # f = open(input_file)
    # global cursor, connection, options
    # options,subtasks=get_options(f)
    # f.close()
    # # cursor,connection = connect_mysql(host=options[0],user=options[1],pw=options[2])
    # tasks = options[4]
    
    # if 'setup' in tasks:
    #     setup(cursor,options[3])
    #     if 'search' in tasks and len(tasks) == 2:
    #         tmp_str = 'Error: Database is empty.'
    #         print_le(tmp_str,'.','Insert values into the database and then search. Aborting due to empty database')
    
    # if 'insert' in tasks:
    #     insert(connection, cursor, options[3],subtasks['insert'])
    
    # if 'search' in tasks:
    #     app.run(debug=True)