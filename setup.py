from setuptools import setup, find_packages

with open("README.md") as readme_file:
    readme = readme_file.read()

setup(
    name='pide',
    version='0.3.1',
    description='A library for petrophysical interpretations of geophysical models.',
    author='Sinan Ozaydin, Lu Li, Utpal Singh, Patrice F. Rey, Maria Constanza Manassero',
    author_email='sinan.ozaydin@protonmail.com',
    url='https://github.com/sinanozaydin/pide',
    python_requires=">3.8",
    license='GNU Lesser General Public License v3 (LGPLv3)',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'matplotlib',
        'scipy',
        'santex',
        'h5py',
        'harmonica',
        'pyproj'
    ],
    keywords = ['petrophysics', 'geodynamic modelling', 'magnetotelluric', 'electrical conductivity', 'seismic velocity']
)
