from setuptools import setup

setup(
    name="m-party",
    version="1.0.0",
    py_modules=["m-party"],
    install_requires=[],
    entry_points={
        "console_scripts": [
            "m-party=m-party:main",
        ],
    },
)