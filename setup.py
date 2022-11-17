from setuptools import find_packages, setup

__version__ = '0.2.0'
URL = 'https://proteinshake.readthedocs.io/en/latest/index.html'

install_requires = [
                    'numpy',
                    'scipy',
                    'biopandas',
                    'scikit-learn',
                    'tqdm',
                    'wget',
                    'requests',
                    'joblib',
                    'rdkit',
                    'tabulate',
                    'fastavro'
                    ]
test_requires = [
    'pytest',
]

setup(
    name='proteinshake',
    version=__version__,
    description='Deep learning ready datasets of 3D protein structures.',
    author = "Tim Kucera, Carlos Oliver, Leslie O'Bray, Dexiong Chen, Karsten Borgwardt",
    author_email = "tim.kucera@bsse.ethz.ch",
    url=URL,
    keywords = ['bioinformatics',
                'deep-learning',
                'pytorch',
                'torch-geometric',
                'computational-biology',
                'macromolecular-structure'],
    python_requires='>=3.7',
    install_requires=install_requires,
    packages=find_packages(),
    include_package_data=True,
)
