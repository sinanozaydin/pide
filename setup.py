from setuptools import setup, find_packages

with open("README.md") as readme_file:
    readme = readme_file.read()

setup(
    name='pide',
    version='0.3.6',
    package_data = {
        'pide' :[
            'pide_src/*csv',
            'pide_src/*json',
            'pide_src/water_sol/*csv',
            'pide_src/water_partitioning/*csv',
            'pide_src/cond_models/*csv',
            'pide_src/cond_models/minerals/*csv',
            'pide_src/cond_models/rocks/*csv'
        ],
    },
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
        'pyproj',
        'netCDF4'
    ],
    keywords = ['petrophysics', 'geodynamic modelling', 'magnetotelluric', 'electrical conductivity', 'seismic velocity']
)
