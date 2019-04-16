from distutils.core import setup
from stk import __version__

setup(name='stk',
      author='Lukas Turcani',
      author_email='lukasturcani93@gmail.com',
      url='https://www.github.com/lukasturcani/stk',
      version=__version__,
      packages=['stk',
                'stk.utilities',
                'stk.molecular',
                'stk.molecular.topologies',
                'stk.molecular.topologies.cage',
                'stk.optimization'],
      install_requires=['networkx',
                        'scipy',
                        'matplotlib',
                        'scikit-learn',
                        'psutil',
                        'pywindowx',
                        'pandas',
                        'seaborn'],
      python_requires='>=3.6')
