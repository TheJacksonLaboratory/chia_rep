from distutils.core import setup
from setuptools.extension import Extension
from Cython.Build import cythonize

extensions = [
    Extension(name="chiapet_rep_util", sources=["chiapet_rep_util.pyx"])
]

setup(
    ext_modules=cythonize(extensions, language_level="3", annotate=True)
)
