from setuptools import setup

setup(
    name="m-party",
    version="0.3.0",
    py_modules=["m_party"],
    install_requires=[
                    "pyyaml", 
                    "hmmer",
                    "pandas", 
                    "openpyxl", 
                    "snakemake", 
                    "clint", 
                    "biopython", 
                    "diamond", 
                    "muscle", 
                    "beautifulsoup4", 
                    "kma", 
                    "tqdm", ],
    entry_points={
        "console_scripts": [
            "m-party=m_party:main",
        ],
    },
)