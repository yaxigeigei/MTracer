from setuptools import setup, find_packages

setup(
    name="msessionexplorer",
    version="0.1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "numpy",
        "pandas",
        "dill",
        "ipython",
    ],
    python_requires=">=3.6",
    author="Duo Xu",
    description="A Python package for managing session data with tables and metadata",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yaxigeigei/MSessionExplorer",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
) 