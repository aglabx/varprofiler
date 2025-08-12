from setuptools import setup, find_packages
import os

def read_requirements():
    """Read requirements from requirements.txt file."""
    req_file = os.path.join(os.path.dirname(__file__), 'requirements.txt')
    if os.path.exists(req_file):
        with open(req_file) as f:
            return [line.strip() for line in f if line.strip() and not line.startswith('#')]
    return ['pandas>=1.3.0', 'matplotlib>=3.3.0']

def read_long_description():
    """Read long description from README.md file."""
    readme_file = os.path.join(os.path.dirname(__file__), 'README.md')
    if os.path.exists(readme_file):
        with open(readme_file, encoding='utf-8') as f:
            return f.read()
    return ''

setup(
    name='varprofiler',
    version='1.0.0',
    author='Your Name',
    author_email='your.email@example.com',
    description='A high-performance tool for analyzing k-mer variability across genomes',
    long_description=read_long_description(),
    long_description_content_type='text/markdown',
    url='https://github.com/yourusername/varprofiler',
    packages=find_packages(),
    py_modules=['plot_chromosomes'],
    install_requires=read_requirements(),
    python_requires='>=3.7',
    entry_points={
        'console_scripts': [
            'varprofiler-plot=plot_chromosomes:main',
        ],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: C++',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS :: MacOS X',
    ],
    keywords='bioinformatics genomics k-mer analysis variation profiling',
    project_urls={
        'Bug Reports': 'https://github.com/yourusername/varprofiler/issues',
        'Source': 'https://github.com/yourusername/varprofiler',
    },
)