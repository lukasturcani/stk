from distutils.core import setup

setup(name='stk',
      author='Lukas Turcani',
      author_email='lukasturcani93@gmail.com',
      url='www.github.com/lukasturcani/stk',
      packages=['stk',
                'stk.convenience_tools',
                'stk.molecular',
                'stk.molecular.topologies',
                'stk.molecular.topologies.cage',
                'stk.optimization'])
