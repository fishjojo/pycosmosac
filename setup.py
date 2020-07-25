#!/usr/bin/env python
from setuptools import setup, find_packages

CLASSIFIERS=[
    'Development Status :: 4 - Beta',
    'Environment :: Console',
    'Environment :: Web Environment',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: MIT',
    'Operating System :: MacOS :: MacOS X',
    'Operating System :: POSIX',
    'Operating System :: Unix',
    'Programming Language :: Python',
    'Topic :: Scientific/Engineering',
]

setup(name='pycosmosac',
      version='0.0.1',
      description='A general implementation of cosmo-sac models',
      author='Xing Zhang',
      author_email='xzhang8@caltech.edu',
      url='https://github.com/fishjojo/pycosmosac',
      license = 'MIT',
      python_requires='>=3.5',
      platforms = ["Linux", "Mac OS-X", "Unix"],
      classifiers = CLASSIFIERS,
      include_package_data=True,
      packages=find_packages(exclude=['*test*', '*example*','*setup.py']),
      install_requires=['numpy>=1.17.2', 'scipy>=1.3.1', 
                        'pandas>=0.25.1','simplejson>=3.17.0',
                        'bs4>=4.8.0'],
      zip_safe=False,
     )
