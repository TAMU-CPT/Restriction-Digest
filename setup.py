#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

requirements = [
    'beautifulsoup4',
    'pyyaml',
    'biopython',
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='dnadigest',
    version='0.0.1',
    description='Command line tool to run restriction digests of DNA sequences',
    author='Stephen Crosby',
    author_email='stcrosby@gmail.com',
    packages=[
        'dnadigest',
    ],
    package_dir={'dnadigest': 'dnadigest'},
    include_package_data=True,
    install_requires=requirements,
    license="GPLv3",
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
