from setuptools import setup, Extension
import numpy

module = Extension(
    'compute_paths_module',
    sources=['compute_paths.c', 'py_compute_paths.c'],
    include_dirs=[numpy.get_include()],
    libraries=["m", "xml2"],
    extra_compile_args=['-O3'],
)

setup(
    name='compute_paths_module',
    version='1.0',
    description='TODO write description.',  # TODO write description
    ext_modules=[module],
)
