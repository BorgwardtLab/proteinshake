from setuptools import find_packages, setup, Extension

__version__ = '0.1.0'
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

tmalign = Extension('TMalign',
    sources = ['tmalign/TMalign.cpp'],
    language = 'c',
    extra_compile_args = ['-static','-O3','-ffast-math','-lm']
)

setup(
    name='proteinshake',
    version=__version__,
    description='Deep learning ready datasets of 3D protein structures.',
    author = "Tim Kucera, Carlos Oliver, Dexiong Chen, Karsten Borgwardt",
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
    ext_modules=[tmalign]
)
