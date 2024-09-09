import os
import re

from setuptools import find_packages, setup

install_requires = \
    ['dask',
     'netcdf4',
     'numpy',
     'pyproj',
     'scipy',
     'xarray>=0.10.0']

here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'pyremap', '__init__.py')) as f:
    init_file = f.read()

version = re.search(r'{}\s*=\s*[(]([^)]*)[)]'.format('__version_info__'),
                    init_file).group(1).replace(', ', '.')

setup(name='pyremap',
      version=version,
      description='Python remapping tools for climate and earth system models',
      url='https://github.com/xylar/pyremap',
      author='Xylar Asay-Davis',
      author_email='xylarstorm@gmail.com',
      license='BSD',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'License :: OSI Approved :: BSD License',
          'Operating System :: OS Independent',
          'Intended Audience :: Science/Research',
          'Programming Language :: Python',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.9',
          'Programming Language :: Python :: 3.10',
          'Programming Language :: Python :: 3.11',
          'Programming Language :: Python :: 3.12',
          'Topic :: Scientific/Engineering',
      ],
      packages=find_packages(),
      package_data={'pyremap.test': ['test*/*', 'test*/*/*']},
      install_requires=install_requires)
