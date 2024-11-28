from setuptools import setup, find_packages

setup(
    name="LaboThapPy",
    version="0.1.0",
    author="Elise Neven, Basile Chaudoir",
    description="Open-source Python library for the modelling and simulation of thermodynamic systems",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",  # Use Markdown format
    url="https://github.com/PyLaboThap/LaboThapPy",
    packages=find_packages(),
    install_requires=[], # Can specify the required packages here (numpy, matplotlib, ...)
    # classifiers=[                    # Metadata for the package
    #     "Programming Language :: Python :: 3",
    #     "License :: OSI Approved :: MIT License",
    #     "Operating System :: OS Independent",
    # ],
    python_requires='>=3.6',         # Minimum Python version required
)

# https://www.youtube.com/watch?v=Mgp6-ZMEcE