#! /usr/bin/env python3
from setuptools import setup

setup(name='EnvClust',
      version='0.3.0',
      packages=['envclust',],
      description='For clustering amplicon sequences using distribution '
          'patterns across samples to inform OTU generation',
      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Natural Language :: English',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3.6',
          'Topic :: Scientific/Engineering :: Bioinformatics',
          'Topic :: Software Development :: Libraries :: Python Modules'
      ],
      keywords='bioinformatics clustering OTU',
      url='https://github.com/cnthornton/envclust/',
      download_url = 'https://github.com/cnthornton/envclust/archive/v0.3.0.tar.gz',
      author='Christopher Thornton',
      author_email='christopher.n.thornton@gmail.com',
      license='GPLv3',
      include_package_data=True,
      zip_safe=False,
      install_requires=['argparse', 'scipy'],
      entry_points={
          'console_scripts': [
              'envclust = envclust.envcluster:main',
          ]
      }
      )
