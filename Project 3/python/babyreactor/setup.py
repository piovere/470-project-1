#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'Click>=6.0',
    # TODO: put package requirements here
    'numpy>=1.18'
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='babyreactor',
    version='0.1.0',
    description="A project implementing a simple 1-d reactor with a reflector",
    long_description=readme + '\n\n' + history,
    author="J.R. Powers-Luhn",
    author_email='jpowersl@vols.utk.edu',
    url='https://github.com/piovere/babyreactor',
    packages=[
        'babyreactor',
    ],
    package_dir={'babyreactor':
                 'babyreactor'},
    entry_points={
        'console_scripts': [
            'babyreactor=babyreactor.cli:main'
        ]
    },
    include_package_data=True,
    install_requires=requirements,
    license="MIT license",
    zip_safe=False,
    keywords='babyreactor',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
