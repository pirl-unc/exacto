#!/usr/bin/env python3


from os.path import dirname, join
from glob import glob
from setuptools import setup, find_packages


DIR = (dirname(__file__) or '.')

setup_args = {}
setup_args.update(
      name='Exacto',
      version='0.1.3',
      description='Expression Quantification and Identification of Augmented Transcripts of Somatic Variant Origin.',
      author='Jin Seok (Andy) Lee',
      author_email='ajslee@unc.edu',
      packages=find_packages(),
      scripts=glob(join(DIR, 'scripts/*.py'))
)

setup(**setup_args)