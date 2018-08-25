#!/usr/bin/env python
import os
from setuptools import setup, find_packages

# https://blog.ionelmc.ro/2014/05/25/python-packaging/
setup(
    name="orpytal",
    version="0.0.2",
    description="Python package for Orbital Mechanics",
    author="Nick LaFarge",
    author_email="nick.lafarge@gmail.com",
    url="http://nicklafarge.space/orpytal",
    download_url="https://github.com/nicklafarge/OrPytal",
    license="MIT",
    keywords=[
        "aero", "aerospace", "engineering",
        "astrodynamics", "orbits", "kepler", "orbital mechanics"
    ],
    python_requires=">=3.5",
    install_requires=[
        "numpy",
        "matplotlib>=2.0",
        "seaborn",
        "matlab",
        "ipython",
        "untangle",
        "pint"
    ],
    extras_require={

    },
    packages=find_packages('src'),
    package_dir={'': 'src'},
    entry_points={
        'console_scripts': [
            'orpytal = orpytal.cli:main'
        ]
    },
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
    long_description=open('README.md').read(),
    include_package_data=True,
    zip_safe=False,
)
