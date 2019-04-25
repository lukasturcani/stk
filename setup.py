from distutils.core import setup
import re
from os.path import join


def get_version():
    with open(join('stk', '__init__.py'), 'r') as f:
        content = f.read()
    p = re.compile(r'^__version__ = [\'"]([^\'\"]*)[\'"]', re.M)
    return p.search(content).group(1)


setup(name='stk',
      author='Lukas Turcani',
      author_email='lukasturcani93@gmail.com',
      url='https://www.github.com/lukasturcani/stk',
      version=get_version(),
      packages=[
          'stk',
          'stk.utilities',
          'stk.molecular',
          'stk.molecular.topologies',
          'stk.molecular.topologies.cage',
          'stk.calculators',
          'stk.calculators.electronic_property',
          'stk.calculators.energy',
          'stk.calculators.optimization'
      ],
      install_requires=[
          'networkx',
          'scipy',
          'matplotlib',
          'scikit-learn',
          'psutil',
          'pywindowx',
          'pandas',
          'seaborn'
       ],
      python_requires='>=3.6')
