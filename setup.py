from setuptools import find_packages, setup
from setuptools.command.install import install
import subprocess

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

class CustomInstall(install):
    """Custom handler for the 'install' command."""
    def run(self):
        print('##################')
        venv = get_virtualenv_path()
        subprocess.check_call(f'g++ -static -O3 -ffast-math -lm -o {venv}/TMalign TMalign.cpp', cwd='.', shell=True)
        super().run()

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
    cmdclass={'install': CustomInstall}
)
