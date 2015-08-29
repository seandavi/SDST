from setuptools import setup, find_packages
import sys, os, glob

version = '0.9.1'

setup(name='SDST',
      version=version,
      description="Useful routines to perform NGS data analysis tasks",
      long_description="""\
""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='',
      author='Sean Davis',
      author_email='seandavi@gmail.com',
      url='https://github.com/seandavi/seqtools',
      license='MIT',
      scripts=glob.glob('scripts/*'),
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          'PyVCF>=0.6.5',
          'pylev'
          #'python-Levenshtein'
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
