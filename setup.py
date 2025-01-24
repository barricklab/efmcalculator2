from setuptools import setup, find_packages

setup(
    name="efmcalculator",
    version="2.1.0",
    python_requires=[3.12],
    install_requires=[
        "pandas",
        "polars==1.19",
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
    packages=find_packages(where="efmcalculator"),
    package_dir={"": "efmcalculator"},
    package_data={"efmcalculator": ["visualization/assets/*"]},
)
