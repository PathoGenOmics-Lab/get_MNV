from setuptools import setup, find_packages

setup(
    name='get_mnv',
    version='0.1',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'PyVCF', 
        'pandas',
        'biopython',
        'argparse'
    ],
    author='Paula Ruiz-Rodriguez, Mireia Coscolla',
    author_email='paula.ruiz-rodriguez@uv.es, mireia.coscolla@uv.es',
    maintainer='Paula Ruiz-Rodriguez',
    maintainer_email='paula.ruiz-rodriguez@uv.es',
    description='A package to use annotate MNVs in VCF files',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/PathoGenOmics-Lab/get_MNV',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10'
    ],
    license='GPLv3',
    python_requires='>=3.6',
    entry_points={
        'console_scripts': [
            'mtbtyper=mtbtyper.main:main',
        ],
    },
)
