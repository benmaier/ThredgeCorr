from setuptools import setup

setup(name='ThredgeCorr',
      version='0.0.2',
      description="Generate instances of George Cantwell's Edge-correlated network models",
      url='https://www.github.com/benmaier/ThredgeCorr',
      author='Benjamin F. Maier',
      author_email='bfmaier@physik.hu-berlin.de',
      license='MIT',
      packages=['ThredgeCorr'],
      install_requires=[
          'numpy>=1.14',
          'scipy>=0.17',
      ],
      zip_safe=False)
