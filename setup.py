from distutils.core import setup

setup(name='stk',
      author='Lukas Turcani',
      author_email='lukasturcani93@gmail.com',
      url='https://www.github.com/lukasturcani/stk',
      version='2018.04.12.0',
      packages=['stk',
                'stk.convenience_tools',
                'stk.molecular',
                'stk.molecular.topologies',
                'stk.molecular.topologies.cage',
                'stk.optimization'],
      install_requires=['networkx',
                        'scipy',
                        'matplotlib',
                        'sklearn',
                        'psutil',
                        'pywindowx'],
      python_requires='>=3.6')
