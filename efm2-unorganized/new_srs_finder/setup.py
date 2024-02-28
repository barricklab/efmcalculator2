from setuptools import setup, find_packages

setup(
    name='short_sequence_finder',
    version='0.1.0',
    packages=find_packages(), 
    install_requires=[
        'Pandas', 
        'Bio',
        'Time'
    ],
    entry_points={
        'console_scripts': [
            'run_short_sequence_finder=short_sequence_finder.short_seq_finder:main',
        ],
    },
    author='Sreya Das',
    author_email='sreyaldas@gmail.com',
    description='This package detects short repeat sequences (6-15) in sequences ',
    url='https://github.com/coolbears/genomescan/short_sequence_finder',
    classifiers=[
        'Programming Language :: Python :: 3',
    ],
)
