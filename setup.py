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
      package_dir={'': 'src'},
      packages=[
          'stk',
          'stk.utilities',
          'stk.molecular',
          'stk.molecular.molecules',
          'stk.molecular.topology_graphs',
          'stk.molecular.topology_graphs.cage',
          'stk.calculators',
          'stk.calculators.electronic_property',
          'stk.calculators.energy',
          'stk.calculators.optimization',
          'stk.calculators.ga'
      ],
      install_requires=[
          'scipy',
          'matplotlib',
          'psutil',
          'pandas',
          'pathos',
          'seaborn',
          'numpy'
       ],
      python_requires='>=3.6')
