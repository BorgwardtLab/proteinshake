from setuptools import find_packages, setup

__version__ = '0.0.3'
URL = 'https://torch-pdb.readthedocs.io/en/latest/index.html'

install_requires = [
                    'numpy',
                    'scipy',
                    'biopandas',
                    'seaborn',
                    'scikit-learn',
                    'tqdm',
                    'wget',
                    'requests',
                    'joblib',
                    'rdkit',
                    'tabulate'
                    ]
test_requires = [
    'pytest',
]

setup(
    name='torch_pdb',
    version=__version__,
    description='PyTorch datasets built from 3D protein structures.',
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

