import setuptools
from distutils.core import setup
from distutils.extension import Extension

extensions = [
    Extension('chia_rep.chiapet_rep_util',
              ['chia_rep/chiapet_rep_util.c']),
]

NAME = 'chia_rep'
VERSION = '0.0.1'

setuptools.setup(

    name=NAME,

    version=VERSION,

    author="Henry Zhang",

    author_email="henrybzhang.99@gmail.com",

    description="A package for finding reproducibility of chIA-PET data.",

    long_description=open('README.md').read(),

    long_description_content_type="text/markdown",

    url="https://github.com/c0ver/chiapet_rep",

    packages=setuptools.find_packages(),

    install_requires=['numpy>=1.17.0',
                      'scipy>=1.3.1',
                      'prettytable>=0.7.2',
                      'pybedgraph>=0.5.22',
                      'matplotlib>=3.1.1'],

    ext_modules=extensions,

    classifiers=[

        "Programming Language :: Python :: 3",

        "License :: OSI Approved :: MIT License",

        "Operating System :: OS Independent",

    ],

)
