from distutils.core import setup
from setuptools.extension import Extension
from Cython.Build import cythonize

extensions = [
    Extension(name="reproducibility_util",
              sources=["reproducibility_util.pyx"])
]

setup(
    ext_modules=cythonize(extensions, language_level="3", annotate=True)
)
