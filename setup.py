from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()


setup(name='pykrev',
      version='0.1',
      description='van Krevelen style analysis and beyond in python',
      long_description = long_description,
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Topic :: Mass spectrometry :: Chemistry',
      ]
      keywords = ['mass spectrometry','van Krevelen','cheminformatics']
      url='https://github.com/Kzra/pykrev',
      author='Ezra Kitson',
      author_email='ezra.kitson@ed.ac.uk',
      license='MIT',
      packages=setuptools.find_packages(),
      install_requires=[
          'numpy',
          'pandas',
	  'matplotlib',
      'Scipy'
      ],
      python_requires='>=3.6'
      zip_safe=False
     )