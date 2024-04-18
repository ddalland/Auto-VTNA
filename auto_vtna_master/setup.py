from setuptools import setup, find_packages
from os import path

# Meta-data of the package.
NAME = 'auto_vtna'
DESCRIPTION = 'A Python package for efficient and automatic VTNA analysis'
EMAIL = 'dd4518@ic.ac.uk'
AUTHOR = 'Daniel Dalland'
VERSION = 1.0
REQUIRED = [
    'numpy',
    'pandas',
    'matplotlib',
    'polyfit',
    'num2words',
    'scikit-learn',
    'scipy'
]

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    LONG_DESCRIPTION = f.read()

setup(
    name=NAME,
    version=VERSION,
    packages=find_packages(),
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    author=AUTHOR,
    author_email=EMAIL,
    python_requires='>=3.6.0',
    long_description_content_type='text/markdown',
    install_requires=REQUIRED,
    include_package_data=True,
    license='MIT',
)
