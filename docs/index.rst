.. ChemBDDB documentation master file, created by
   sphinx-quickstart on Wed Oct 30 16:45:02 2019.

ChemBDDB: A Big Data Database Toolkit for Chemical and Materials Data Storage
==============================================================================

ChemBDDB is a big data database toolkit which facilitates the efficient management and sharing of chemical and materials data. It can setup a database, populate it with the data provided by the user and enable the user to search the database, visualize each molecule in the database and share the data via the web.

Program Version: 0.1

Release Date: July 5, 2019

With contributions by: 

- Shirish Sivraj (UB): Molecule visualization 
- Supriya Agarwal (UB): Initial framework

Code Design:
++++++++++++

ChemBDDB is developed in the Python 3 programming language and uses MySQL, Flask to set up a SQL database using pymysql. It uses OpenBabel and its Python extension, Pybel for handling molecules. The development follows a strictly modular and object-oriented design to make the overall code as flexible and versatile as possible.

Installation and Dependencies:
++++++++++++++++++++++++++++++

Installing MySQL is a prerequisite for running ChemBDDB. General instructions:

    - For Windows:

        - install mysql server: https://dev.mysql.com/downloads/mysql/ 1.1 create a simple password for the root user and make a note of it
        
        - Add the path to \MySQL Server 8.0\bin to environment/path variables

        - It is highly recommended that a virtual environment is used to run ChemBDDB. The virtual environment and ChemBDDB and its dependencies can be installed as:

            .. code:: bash

                conda create --name my_chembddb_env python=3.7
                source activate my_chembddb_env
                conda install -c openbabel openbabel
                pip install chembddb

        - You can test the installation with:

            .. code:: bash

                pytest -v

Documentation:
++++++++++++++

.. toctree::
    chembddb
    chembddb.python_module
    chembddb.python_module_tutorial

Citation:
+++++++++

- Sonpal, A.; Agrawal, S.; Sivaraj, S.; Hachmann, ChemBDDB – A Big Data Database Toolkit for Chemical and Materials Data Storage. 2019;https://github.com/hachmannlab/chembddb

- J. Hachmann, M.A.F. Afzal, M. Haghighatlari, Y. Pal, Building and Deploying a Cyberinfrastructure for the Data-Driven Design of Chemical Systems and the Exploration of Chemical Space, Mol. Simul. 44 (2018), 921-929. DOI: 10.1080/08927022.2018.1471692

Acknowledgement:
++++++++++++++++

ChemBDDB is based upon work supported by the U.S. National Science Foundation under grant #OAC-1751161. It was also supported by start-up funds provided by UB's School of Engineering and Applied Science and UB's Department of Chemical and Biological Engineering, the New York State Center of Excellence in Materials Informatics through seed grant #1140384-8-75163, and the U.S. Department of Energy under grant #DE-SC0017193.

License and Copyright
+++++++++++++++++++++

- ChemBDDB is distributed under 3-Clause BSD License (https://opensource.org/licenses/BSD-3-Clause).
- ChemBDDB is copyright (C) 2015-2018 Johannes Hachmann and Aditya Sonpal, all rights reserved. 
- (C) 2015-2018 Johannes Hachmann, Aditya Sonpal
- University at Buffalo - The State University of New York (UB)

Contact: 
- hachmann@buffalo.edu
- adityaso@buffalo.edu
- http://hachmannlab.cbe.buffalo.edu