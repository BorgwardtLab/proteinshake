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

def compile_and_install_software():
    """Used the subprocess module to compile/install the C software."""
    src_path = './some_c_package/'

    # compile the software
    cmd = "./configure CFLAGS='-03 -w -fPIC'"
    venv = get_virtualenv_path()
    if venv:
        cmd += ' --prefix=' + os.path.abspath(venv)
    subprocess.check_call(cmd, cwd=src_path, shell=True)

    # install the software (into the virtualenv bin dir if present)
    subprocess.check_call('make install', cwd=src_path, shell=True)


class CustomInstall(install):
    """Custom handler for the 'install' command."""
    def run(self):
        subprocess.check_call('g++ -static -O3 -ffast-math -lm -o TMalign TMalign.cpp', cwd='.', shell=True)
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
