from setuptools import setup, find_packages

setup(
    name="efmcalculator",
    version="2.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "polars",
        "progress",
        "biopython",
        "statsmodels",
        "rich",
        "pyarrow",
        "bokeh==2.4.3",
        "numpy<2",
        "streamlit-extras",
    ],
    entry_points={
        "console_scripts": [
            "efmcalculator-webapp=efmcalculator.webapp.bootstrap_streamlit:bootstrap_streamlit",
        ],
    },
    author="placeholder",
    author_email="placeholder",
    description="placeholder",
    url="https://github.com/barricklab/efm-calculator2",
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    include_package_data=True,
)
