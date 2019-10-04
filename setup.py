from setuptools import setup, find_packages

version = '0.0.1'

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
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Topic :: Scientific/Engineering',
      ],
      packages=find_packages(),
      package_data={'pyremap.test': ['test*/*', 'test*/*/*']})
