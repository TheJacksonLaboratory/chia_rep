from distutils.core import setup
from setuptools.extension import Extension
from Cython.Build import cythonize

extensions = [
    Extension(name="chia_rep_util", sources=["chia_rep_util.pyx"])
]

setup(
    ext_modules=cythonize(extensions, language_level="3", annotate=True)
)
