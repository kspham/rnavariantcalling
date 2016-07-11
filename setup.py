from setuptools import setup

setup(name='rnavariantcalling',
      version='0.1',
      description='Efficiency tool for calling SNP variant of RNA sequence',
      url='http://github.com/kspham/rnavariantcalling',
      author='Ptdtan',
      author_email='ptdtan@gmail.com',
      license='MIT',
      packages=['rnavariantcalling'],

      install_requires=['pyyaml'],
      scripts=['src/multithread.py', 'src/rnavariantcalling.py', 'scripts/filter', 'scripts/filter2', 'scripts/filter3',
	  			'scripts/vcffirstheader', 'scripts/vcfstreamsort', 'scripts/vcfuniq', 'scripts/firstheader', 'scripts/fby', 'src/freebayes_pool.py'],
      zip_safe=False)
