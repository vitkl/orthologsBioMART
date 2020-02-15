from setuptools import setup
from setuptools import find_packages

import sys

def setup_package():
  install_requires = ['pybiomart', 'pandas', 'tqdm']
  metadata = dict(
      name = 'pyorthomap',
      version = '0.1',
      description = 'pyorthomap: Map orthologous genes using ENSEMBL BioMart and pybiomart package',
      url = 'https://github.com/vitkl/orthologsBioMART',
      author = 'Vitalii Kleshchevnikov',
      author_email = 'vitalii.kleshchevnikov@sanger.ac.uk',
      license = 'Apache License, Version 2.0',
      packages = find_packages(),
      install_requires = install_requires
    )

  setup(**metadata)

if __name__ == '__main__':
  if sys.version_info < (2,7):
    sys.exit('Sorry, Python < 2.7 is not supported')

  setup_package()
