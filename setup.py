"""
openff-fragmenter

Fragment molecules for quantum mechanics torsion scans.
"""

import sys

from setuptools import find_namespace_packages, setup

import versioneer

short_description = __doc__.split("\n")

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {"pytest", "test", "ptr"}.intersection(sys.argv)
pytest_runner = ["pytest-runner"] if needs_pytest else []

try:
    with open("README.md") as handle:
        long_description = handle.read()
except OSError:
    long_description = "\n".join(short_description[2:])


setup(
    name="openff-fragmenter",
    author="Chaya D. Stern",
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="MIT",
    packages=find_namespace_packages(include=["openff.*"]),
    include_package_data=True,
    setup_requires=[] + pytest_runner,
    python_requires=">=3.6",
)
