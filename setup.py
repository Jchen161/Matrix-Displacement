from setuptools import setup, find_packages

setup(
    name="matrix-displacement",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.21.0",
        "matplotlib>=3.4.0",
    ],
    author="Jinxin Chen",
    author_email="jchen161@stevens.edu",
    description="Matrix Displacement Method implementation for structural analysis",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/JinxinDev/Matrix-Displacement.git",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Education",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Education",
    ],
    python_requires=">=3.8",
)