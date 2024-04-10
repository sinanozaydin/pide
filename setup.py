from setuptools import setup, find_packages

setup(
    name='pide',
    version='0.1.1',
    description='A simple Python library',
    author='Sinan Ozaydin',
    author_email='sinan.ozaydin@protonmail.com',
    url='https://github.com/sinanozaydin/pide',
    license='MIT',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'matplotlib',
        'scipy',
        'satex',
        'h5py'
    ],
)