#!/usr/bin/env python

from pathlib import Path
from setuptools import setup, find_packages

# Read the README for the long description if available
here = Path(__file__).parent
long_description = ""
readme = here / "README.md"
if readme.exists():
    long_description = readme.read_text(encoding="utf-8")

# Top-level modules (single-file scripts) that live in the project root
# These are included so users can import them after installation.
py_modules = [
    "Main",
    "start",
    "Reciprocal",
    "Report",
    "dom",
    "initblast",
]

setup(
    name="pevolution",
    version="1.1.0",
    description=(
        "Pipeline to find, align, and build trees for putatively related proteins."
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Mihir Sarwade",
    author_email="mihir.sarwade@gmail.com",
    packages=find_packages(exclude=["test", "venv1"]),
    py_modules=py_modules,
    include_package_data=True,
    package_data={"": ["Data/*", "Data/**/*"]},
    install_requires=["psutil>=5.6.3", "biopython>=1.86"],
    python_requires=">=3.10",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    # Intentionally not creating a console_scripts entry point because
    # existing top-level modules expose `main(argv)` (not a zero-arg callable).
    # If you'd like a CLI entry point, I can add a small wrapper module
    # (for example `cli.py`) which invokes `Main.main(sys.argv[1:])` and
    # then add `entry_points={'console_scripts': ['pevolution=cli:main']}`.
)
