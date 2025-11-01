from setuptools import setup, find_packages

setup(
    name="mosaicprot",
    version="0.1.13",
    packages=find_packages(),
    install_requires=[
        # List your dependencies here
    "biopython==1.85",
    "pandas==2.2.3"
    ],
    entry_points={
        "console_scripts": [
            "mosaicprot=mosaicprot.cli:run",
        ],
    },
    author="Umut Cakir & Ali Yurtseven",
    author_email="hpumut@gmail.com",
    description="MosaicProt CLI: Detect ORFs, Find Alternative ORFs, Generate Mosaic Proteins",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url= "https://github.com/umutcakir/MosaicProt",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
