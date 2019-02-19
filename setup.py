from setuptools import setup

setup(name='ThredgeCorr',
      version='0.0.4',
      description="Generate instances of the thresholded locally-correlated edge weights network model.",
      url='https://www.github.com/benmaier/ThredgeCorr',
      author='Benjamin F. Maier',
      author_email='bfmaier@physik.hu-berlin.de',
      license='MIT',
      packages=['ThredgeCorr'],
      install_requires=[
          'numpy>=1.14',
          'scipy>=0.17',
          'bfmplot>=0.0.4',
          'mpmath',

      ],
      zip_safe=False)
