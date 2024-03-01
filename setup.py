from setuptools import setup, find_packages

setup(
    name='Evolutionary Failure Mode Calculator',
    version='2.1.0',
    packages=find_packages(), 
    install_requires=[
        'pandas', 
        'biopython'
    ],
    entry_points={
        'console_scripts': [
            'run_short_sequence_finder=short_sequence_finder.short_seq_finder:main',
        ],
    },
    author='placeholder',
    author_email='placeholder',
    description='placeholder',
    url='https://github.com/barricklab/efm-calculator2',
    classifiers=[
        'Programming Language :: Python :: 3',
    ],
)
