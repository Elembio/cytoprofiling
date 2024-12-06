from setuptools import setup, find_packages

setup(
    name="cytoprofiling",
    version="1.0.0",
    description="Package for processing data from Element Biosciences Teton cytoprofiling platform",
    author="Ryan Kelley",
    author_email="ryan.kelley@elembio.com",
    packages=find_packages(exclude=("tests", "docs", "examples")),
    install_requires=[
        'numpy',
        'pandas',
        'typing'
    ],
    extras_require = {
        "scanpy": ["scanpy", "anndata"]
    }
)
