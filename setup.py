from setuptools import find_packages, setup

setup(
    name="md",
    packages=find_packages(),
    install_requires=["docopt", "packaging", "PeptideBuilder"],
)
