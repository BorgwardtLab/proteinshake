from setuptools import find_packages, setup
from setuptools.command.install import install
from setuptools.command.develop import develop
from setuptools.command.egg_info import egg_info
import subprocess, sys

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

def get_virtualenv_path():
    """Used to work out path to install compiled binaries to."""
    if hasattr(sys, 'real_prefix'):
        return sys.prefix
    if hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix:
        return sys.prefix
    if 'conda' in sys.prefix:
        return sys.prefix
    return None

def install_tm():
    venv = get_virtualenv_path()
    print('#',venv)
    subprocess.check_call(f'g++ -static -O3 -ffast-math -lm -o {venv}/TMalign TMalign.cpp', cwd='.', shell=True)


class CustomInstallCommand(install):
    def run(self):
        install_tm()
        install.run(self)

class CustomDevelopCommand(develop):
    def run(self):
        install_tm()
        develop.run(self)

class CustomEggInfoCommand(egg_info):
    def run(self):
        install_tm()
        egg_info.run(self)


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
    cmdclass={
        'install': CustomInstallCommand,
        'develop': CustomDevelopCommand,
        'egg_info': CustomEggInfoCommand,
    }
)
