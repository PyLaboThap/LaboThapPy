# from setuptools import setup, find_packages

# setup(
#     name='LaboThapPy',          # Nom du package
#     version='1.0',              # Version du package
#     packages=find_packages()    # Inclut tous les sous-packages automatiquement
# )

# https://www.youtube.com/watch?v=Mgp6-ZMEcE
from setuptools import setup, find_packages

setup(
    name="LaboThapPy",  # Library's name
    version="0.1.0",           # Initial version
    author="Elise Neven",
    author_email="elise.neven@uliege.be",
    description="A library to modelize and simulate thermodynamic systems.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown", # README file type
    url="https://github.com/PyLaboThap/LaboThapPy",  # GitHub repo URL
    packages=find_packages(),  # Automatically find packages in your_library/
    install_requires=[         # Dependencies your library needs
        "numpy>=2.1.1",
        "pandas>=1.3.0",
        "CoolProp==6.6.0",
        "pandas==2.2.3",
        "scipy==1.14.1"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
