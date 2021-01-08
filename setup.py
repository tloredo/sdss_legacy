import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
  long_description = fh.read()

setuptools.setup(
    name="sdss_legacy",  # Replace with your own username
    version="2021.1a1",
    author="Tom Loredo",
    author_email="loredo@astro.cornell.edu",
    description="Access SDSS-I/SDSS-II MGS and QSO catalog and spectral data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="TBD",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD 3-clause License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
)
