import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ot_tractability_pipeline",
    version="0.0.1",
    author="Chris Radoux",
    author_email="cradoux@ebi.ac.uk",
    description="Python script for assigning genes to tractability buckets",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Cradoux/tractability_pipeline",
    packages=setuptools.find_packages(),
    install_requires=[i.strip() for i in open("requirements.txt").readlines()],
    include_package_data=True,
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)