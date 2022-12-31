import setuptools

#with open("README.md", "r") as fh:
#    long_description = fh.read()

setuptools.setup(
    name="pymtp",
    version="0.1",
    author="Hao Gao, Junjie Wang",
    author_email="gaaooh@126.com",
    description="A simple Python package for MTP calculator.",
    #long_description=long_description,
    #long_description_content_type="text/markdown",
    #url="https://github.com/gooaah/pycqg",
    # include_package_data=True,
    # exclude_package_date={'':['.gitignore']},
    packages=setuptools.find_packages(),
    install_requires=[
        'numpy',
        'ase',
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)
