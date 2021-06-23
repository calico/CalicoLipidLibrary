from setuptools import setup, find_packages

VERSION = '0.0.3'
DESCRIPTION = 'In silico lipid fragmentation generator (MS/MS spectra)'
LONG_DESCRIPTION = 'In silico lipid fragmentation generator (MS/MS spectra), based on chemical standards and scientific literature.'

# Setting up
setup(
    name="calicolipidlibrary",
    version=VERSION,
    author="Bryson Bennett",
    author_email="bryson@calicolabs.com",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    install_requires=[],
    keywords=['python', 'calico', 'lipid', 'MS/MS', 'in silico', 'fragmentation', 'spectral library', 'spectrum', 'MS2'],
    classifiers=[
        "Environment :: Console",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 2",
    ]
)