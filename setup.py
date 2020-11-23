import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pykrev", 
    version="0.0.1",
    author="Ezra Kitson",
    author_email="ezrakitson@ed.ac.uk",
    description="van Krevelen analysis and beyond in Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kzra/pykrev",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)