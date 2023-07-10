from setuptools import setup

setup(name='ThredgeCorr',
      version='0.1.0',
      description="Generate instances of the thresholded locally-correlated edge weights network model.",
      url='https://www.github.com/benmaier/ThredgeCorr',
      author='Benjamin F. Maier, George C. Cantwell, Guillaume St-Onge',
      author_email='bfmaier@physik.hu-berlin.de',
      license='MIT',
      packages=['ThredgeCorr'],
      install_requires=[
          'numpy>=1.14',
          'scipy>=0.17',
          'mpmath>=1.3',
          'networkx>=2.0',
      ],
      zip_safe=False)
