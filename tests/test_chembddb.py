"""
Unit and regression test for chembddb.py
"""

from chembddb import chembddb
import pytest
import pymysql
import sys

def test_chembddb_imported():
    assert 'chembddb' in sys.modules

@pytest.mark.parametrize("host,user,pw,res",[('127.0.0.1','root','',pymysql.cursors.Cursor),('127.0.0.1','root','wrongpw',str)])
def test_connect_mysql(host,user,pw,res):
    """Testing mysql connect function"""
    cur,all_dbs=chembddb.connect_mysql(host,user,pw)
    assert type(cur) is res

@pytest.mark.parametrize("host,user,pw,db,res",[('127.0.0.1','root','','ben','Success'),('127.0.0.1','root','','ben','Failed! Database already exists.')])
def test_setup(host,user,pw,db,res):
    result = chembddb.setup(host=host,user=user,pw=pw,db=db)
    assert result == res

@pytest.mark.parametrize("host,user,pw,db,smi_col,mol_identifier,conf_file,data_file,res",[('127.0.0.1','root','','ben','smiles','name','tests/test-files/config_hsp.csv','tests/test-files/benzene.csv','Successfully entered the data into the database'),('127.0.0.1','root','','newdb','smiles','name','tests/test-files/config_hsp.csv','tests/test-files/benzene.csv','Failed! Database does not exists.'),('127.0.0.1','root','','ben','smiles','name','tests/test-files/config_hsp.csv','tests/test-files/benzene.csv','Failed! Duplicate entries for all molecules exist.'),('127.0.0.1','root','','ben','smiles','name','tests/test-files/config_hsp.csv','tests/test-files/r2_length7_properties.csv','Failed! Property(s) in config file do not exist in data.'),('127.0.0.1','root','','ben','smiles','noname','tests/test-files/config_hsp.csv','tests/test-files/benzene.csv','Failed! Identifier(s) listed do not exist in data.'),('127.0.0.1','root','','ben','smiles','common_name','tests/test-files/config_hsp.csv','tests/test-files/duplicate_free_hansen_q.csv','A few molecules were not entered due to duplicate entries but the database was successfully populated with the rest'),('127.0.0.1','root','','ben','smiles','common_name','tests/test-files/configuration_hsp.csv','tests/test-files/duplicate_free_hansen_q.csv','Failed! Config file does not exist in the path specified'),('127.0.0.1','root','','ben','smiles','common_name','tests/test-files/config_hsp.csv','tests/test-files/duplicate.csv','Failed! Data file does not exist in the path specified.')])
def test_insert(host, user, pw, db, smi_col, mol_identifier, conf_file, data_file, res):
    result = chembddb.insert(host=host, user=user, pw=pw, db=db, smi_col=smi_col, mol_identifier=mol_identifier, conf_file=conf_file, data_file=data_file)
    assert result == res

@pytest.mark.parametrize("host,user,pw,db,res",[('127.0.0.1','root','','ben','Successfully deleted the database'),('127.0.0.1','root','','ben','Failed! database does not exist')])
def test_delete(host, user, pw, db, res):
    result = chembddb.delete(host=host, user=user, pw=pw, db=db)
    assert result == res

