"""
Unit and regression test for chembddb.py
"""

from chembddb import chembddb
import pytest
import pymysql

@pytest.mark.parametrize("host,user,pw,res",[('127.0.0.1','root','aditya123',pymysql.cursors.Cursor),('127.0.0.1','root','wrongpw',str)])
def test_connect_mysql(host,user,pw,res):
    """Testing mysql connect function"""
    cur,all_dbs=chembddb.connect_mysql(host,user,pw)
    assert type(cur) is res

@pytest.mark.parametrize("host,user,pw,name,res",[('127.0.0.1','root','aditya123','nn1',True),('127.0.0.1','root','aditya123','nn1',False)])
def test_setup(host,user,pw,name,res):
    result = chembddb.setup(host=host,user=user,pw=pw,name=name)
    assert result == res