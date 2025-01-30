from setuptools import setup

setup(
    name="m-party",
    version="1.0.0",
    py_modules=["m-party"],
    install_requires=["upimapi", 
                    "pyyaml", 
                    "hmmer", 
                    "t-coffee", 
                    "cd-hit", 
                    "pandas", 
                    "openpyxl", 
                    "snakemake", 
                    "clint", 
                    "blast", 
                    "diamond", 
                    "muscle", 
                    "beautifulsoup4", 
                    "kma", 
                    "tqdm"],
    entry_points={
        "console_scripts": [
            "m-party=m-party:main",
        ],
    },
)