import subprocesses
import pymysql 
def check_db(cur, con, db):
    

def backup(user, pw, db, filename):
    """Set docstring here.

    Parameters
    ----------
    user: str
        username for mysql server
    pw: str
        password for mysql server
    db: str
        name of the database that needs to be backed up
    filename: str
        full path and name of the backup file

    """ 
    try:
        s= 'mysqldump -h localhost -P 3306 -u '+user+' -p'+pw+ ' '+db+' --single-transaction  > ' + filename
        subprocesses.Popen(s,shell=True)
        return 'done'
    except:
        return 'error'
