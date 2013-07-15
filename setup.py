from setuptools import setup, Extension, find_packages
import sys

desc = 'Tools for phasing and separating SMRT sequencing data'

if ("install" in sys.argv) and sys.version_info < (2, 7, 0):
    raise SystemExit("PhasingTools requires Python 2.7")

setup(
    name = 'PhasingTools',
    version='0.1.0',
    author='Brett Bowman',
    author_email='bbowman@pacificbiosciences.com',
    url='https://github.com/bnbowman/PhasingTools',
    description=desc,
    license=open('LICENSES.txt').read(),
    packages = find_packages('src'),
    package_dir = {'':'src'},
    zip_safe = False,
    install_requires=[]
    )
