import setuptools
from distutils.extension import Extension
from Cython.Build import cythonize

extensions = [
    Extension('chia_rep.chia_rep_util',
              ['chia_rep/chia_rep_util.pyx']),
]

NAME = 'ChIA-Rep'
VERSION = '1.0.0'

setuptools.setup(

    name=NAME,

    version=VERSION,

    author="Henry Zhang",

    author_email="henrybzhang.99@gmail.com",

    description="A package for measuring reproducibility of ChIA-PET data.",

    long_description=open('README.md').read(),

    long_description_content_type="text/markdown",

    url="https://github.com/c0ver/chia_rep",

    packages=['chia_rep'],

    include_package_data=True,

    install_requires=['numpy>=1.17.0',
                      'scipy>=1.3.1',
                      'prettytable>=0.7.2',
                      'pybedgraph>=0.5.40',
                      'matplotlib>=3.1.1'],

    ext_modules=cythonize(extensions, language_level=3),

    license="MIT",

    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],

)
