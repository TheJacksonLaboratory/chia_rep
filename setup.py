import setuptools
from distutils.core import setup
from distutils.extension import Extension

extensions = [
    Extension('chiapet_rep.chiapet_rep_util',
              ['chiapet_rep/chiapet_rep_util.c']),
]

NAME = 'chiapet_rep'
VERSION = '0.0.01'

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

    install_requires=['numpy', 'scipy', 'prettytable', 'pybedgraph', 'matplotlib'],

    ext_modules=extensions,

    classifiers=[

        "Programming Language :: Python :: 3",

        "License :: OSI Approved :: MIT License",

        "Operating System :: OS Independent",

    ],

)
