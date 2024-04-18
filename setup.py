from setuptools import setup, find_packages

setup(
    name='pide',
    version='0.1.1',
    description='A library for petrophysical interpretation.',
    author='Sinan Ozaydin',
    author_email='sinan.ozaydin@protonmail.com',
    url='https://github.com/sinanozaydin/pide',
    license='MIT',
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
)
