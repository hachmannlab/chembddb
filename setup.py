import setuptools
from os import path
import chembddb

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md')) as f:
    long_description = f.read()

if __name__ == "__main__":
    setuptools.setup(
        name='chembddb',
        version=chembddb.__version__,
        author='Aditya Sonpal, Johannes Hachmann',
        author_email='adityaso@buffalo.edu, hachmann@buffalo.edu',
        # url='https://github.com/hachmannlab/chemml',
        project_urls={
            'Source': 'https://github.com/hachmannlab/chembddb',
            'url': 'https://hachmannlab.github.io/chembddb/'
        },
        description=
        'chembddb is a big data database generator.',
        long_description=long_description,
        # scripts=['lib/chembddbshell'],
        include_package_data=True,
        keywords=[
            'Big data', 'SQL database',
            'Web-framework','Materials Science', 
        ],
        license='BSD-3C',
        packages=setuptools.find_packages(),
        scripts = ['lib/chembddbshell'],
        install_requires=[
            'numpy', 'pandas',
            'flask','flask-mysql',
            'pymysql',
        ],
        extras_require={
            'docs': [
                'sphinx',
                'sphinxcontrib-napoleon',
                'sphinx_rtd_theme',
                'numpydoc',
            ],
            'tests': [
                'pytest',
                'pytest-cov',
                'pytest-pep8',
                'tox',
            ],
        },
        tests_require=[
            'pytest',
            'pytest-cov',
            'pytest-pep8',
            'tox',
        ],
        classifiers=[
            'Development Status :: 4 - Beta',
            'Natural Language :: English',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.6',
        ],
        zip_safe=False,
    )
