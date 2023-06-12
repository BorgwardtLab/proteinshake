import datetime
import proteinshake
from sphinxawesome_theme import ThemeOptions

author = "Tim Kucera, Carlos Oliver, Dexiong Chen, Karsten Borgwardt"
project = 'ProteinShake'
version = "0.0.1"
copyright = f'{datetime.datetime.now().year}, {author}'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'myst_nb',
    'sphinxawesome_theme.highlighting'
]
html_permalinks = False
html_static_path = ['_static']
html_css_files = [
    'style.css',
]

html_theme = "sphinxawesome_theme"
#html_logo = ('https://raw.githubusercontent.com/BorgwardtLab/proteinshake/main/docs/images/logo_docs.png')
html_favicon= ('https://raw.githubusercontent.com/BorgwardtLab/proteinshake/main/docs/images/favicon.ico')

add_module_names = False
autodoc_member_order = 'bysource'

intersphinx_mapping = {
    'python': ('https://docs.python.org/', None),
    'numpy': ('http://docs.scipy.org/doc/numpy', None),
    'pandas': ('http://pandas.pydata.org/pandas-docs/dev', None),
    'torch': ('https://pytorch.org/docs/master', None),
}

nb_execution_timeout = 60*10


def setup(app):
    def rst_jinja_render(app, _, source):
        rst_context = {'proteinshake': proteinshake}
        source[0] = app.builder.templates.render_string(source[0], rst_context)

    app.connect('source-read', rst_jinja_render)

