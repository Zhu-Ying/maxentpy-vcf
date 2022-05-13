import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="maxentpy-vcf",
    version="v0.0.1",
    author="Ying Zhu",
    author_email="win19890412@163.com",
    description="A tool for MaxEntScan base on VCF",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'pyfaidx', 'pandas', 'PyVCF'
    ],
    package_data={
        "maxentpy_vcf": ['maxentpy/data/*'],
    },
    scripts=['maxentpy-vcf.py']
)
