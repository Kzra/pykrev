from setuptools import setup

setup(name='pykrev',
      version='0.1',
      description='van Krevelen style analysis and beyond in python',
      long_description = 'Please visit the GitHub repository for documentation and a user guide',
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Topic :: Mass spectrometry :: Chemistry',
      ]
      keywords = ['mass spectrometry','van Krevelen','cheminformatics']
      url='https://github.com/Kzra/PyKrev',
      author='Ezra Kitson',
      author_email='ezra.kitson@ed.ac.uk',
      license='MIT',
      packages=['pykrev'],
      install_requires=[
          'numpy',
          'pandas',
	  'matplotlib',
      'Scipy'
      ],
      zip_safe=False
     )