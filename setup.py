import setuptools

setuptools.setup(
    name="pelevs",
    version="0.0.1",
    author="Ignasi Puch-Giner",
    author_email="ignasi.puch.giner@gmail.com",
    description="A Python package to automatize large dataset VS campaign.",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)