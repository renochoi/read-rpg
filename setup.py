## setup.py

from glob import glob
from os.path import basename, splitext
from setuptools import find_packages, setup

setup(
    name='rpg-tools',
    version='0.1',
    author='Reno Choi',
    author_email='renochoi@gmail.com',
    #long_description=read('README.md'),
    python_requires='>=3.6',
    #install_requires=['numpy'],
    package_data=setuptools.find_packages()
    #dependency_links = [], ## 최신 패키지를 설치하는 경우 
    description='my project',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
)