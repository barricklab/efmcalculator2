from setuptools import setup, find_packages

setup(
    name="efmcalculator",
    version="2.1.0",
    packages=find_packages(),
    install_requires=["pandas", "progress", "biopython", "statsmodels", "rich", "pyarrow"],
    entry_points={
        "console_scripts": [
            "efmcalculator=efmcalculator.efmcalculator:_main",
        ],
    },
    author="placeholder",
    author_email="placeholder",
    description="placeholder",
    url="https://github.com/barricklab/efm-calculator2",
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
)
