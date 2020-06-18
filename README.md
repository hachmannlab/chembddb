[![license](http://img.shields.io/badge/license-BSD-blue.svg?style=flat)](https://github.com/hachmannlab/chemml/blob/master/LICENSE)
[![Build Status](https://travis-ci.org/hachmannlab/chembddb.svg?branch=master)](https://travis-ci.org/hachmannlab/chembddb)
[![Documentation Status](https://readthedocs.org/projects/chembddb/badge/?version=latest)](https://chembddb.readthedocs.io/en/latest/?badge=latest)
# ChemBDDB
ChemBDDB is a big data database toolkit which facilitates the efficient management and sharing of chemical and materials data. It can setup a database, populate it with the data provided by the user and enable the user to search the database, visualize each molecule in the database and share the data via the web.

<p align="center">
  <img align="middle" src="./docs/images/chembddb_logo.png" alt="chembddb" width="400px" class="center">
 </p>

Program Version: 0.1

Release Date: July 5, 2019

With contributions by:<br/> 
Shirish Sivraj (UB): Molecule visualization <br/> 
Supriya Agarwal (UB): Initial framework

## Code Design:
ChemBDDB is developed in the Python 3 programming language and uses MySQL, Flask to set up a SQL database using pymysql. It uses OpenBabel and its Python extension, Pybel for handling molecules. The development follows a strictly modular and object-oriented design to make the overall code as flexible and versatile as possible. 

## Documentation:
ChemBDDB documentation can be found here https://chembddb.readthedocs.io/en/latest/chembddb.html

## Installation and Dependencies:
Installing MySQL is a prerequisite for running ChemBDDB. General instructions:

For Windows
1. install mysql server: https://dev.mysql.com/downloads/mysql/ 
    1.1 create a simple password for the root user and make a note of it
2. Add the path to \MySQL Server 8.0\bin to environment/path variables

It is highly recommended that a virtual environment is used to run ChemBDDB. The virtual environment and ChemBDDB and its dependencies can be installed as:

    conda create --name my_chembddb_env python=3.7
    source activate my_chembddb_env
    conda install -c openbabel openbabel
    pip install chembddb

## Citation:
Please cite the use of ChemBDDB as:

    (1) Sonpal, A.; Agrawal, S.; Sivaraj, S.; Hachmann, J.Chem-BDDB– A Big Data Database Toolkit for Chemical and Materials Data Storage. 2019;https://github.com/hachmannlab/chembddb.
    (2) J. Hachmann, M.A.F. Afzal, M. Haghighatlari, Y. Pal, Building and Deploying a Cyberinfrastructure for the Data-Driven Design of Chemical Systems and the Exploration of Chemical Space, Mol. Simul. 44 (2018), 921-929. DOI: 10.1080/08927022.2018.1471692

## Acknowledgement
ChemBDDB is based upon work supported by the U.S. National Science Foundation under grant #OAC-1751161. It was also supported by start-up funds provided by UB's School of Engineering and Applied Science and UB's Department of Chemical and Biological Engineering, the New York State Center of Excellence in Materials Informatics through seed grant #1140384-8-75163, and the U.S. Department of Energy under grant #DE-SC0017193.

## License and Copyright:
ChemBDDB is copyright (C) 2015-2018 Johannes Hachmann and Aditya Sonpal, all rights reserved. 
ChemBDDB is distributed under 3-Clause BSD License (https://opensource.org/licenses/BSD-3-Clause).

(C) 2015-2018 Johannes Hachmann, Aditya Sonpal
University at Buffalo - The State University of New York (UB)
Contact: hachmann@buffalo.edu, adityaso@buffalo.edu
http://hachmannlab.cbe.buffalo.edu
