import setuptools
from distutils.extension import Extension

USE_CYTHON = False
try:
    from Cython.Distutils import build_ext
    USE_CYTHON = True
except ImportError:
    pass

ext = 'pyx' if USE_CYTHON else 'c'
extensions = [
    Extension('chia_rep.chia_rep_util',
              [f'chia_rep/chia_rep_util.{ext}']),
]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions, language_level=3)

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

    ext_modules=extensions,

    license="MIT",

    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],

)
