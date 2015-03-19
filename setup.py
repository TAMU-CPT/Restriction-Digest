#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

readme = open('README.rst').read()
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

requirements = [
    'beautifulsoup4',
    'pyyaml',
    'biopython',
    'svgwrite',
]

test_requirements = [
    'nose',
]

setup(
    name='dnadigest',
    version='0.1.2',
    description='Library+Tool to run restriction digests of DNA sequences',
    long_description=readme + '\n\n' + history,
    author='Stephen Crosby',
    author_email='stcrosby@gmail.com',
    url='https://github.com/erasche/restriction-digest/',
    packages=[
        'dnadigest',
    ],
    scripts=[
        'bin/digest_dna.py',
        'bin/graphic_drawer.py',
    ],
    package_dir={'dnadigest': 'dnadigest'},
    include_package_data=True,
    install_requires=requirements,
    license="GPLv3",
    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    test_suite='test',
    tests_require=test_requirements
)
