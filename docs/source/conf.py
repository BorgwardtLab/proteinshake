import datetime

import proteinshake

author = "Kucera, Oliver, O'Bray, Chen"
project = 'proteinshake'
version = "0.0.1"
copyright = f'{datetime.datetime.now().year}, {author}'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
]

html_theme = 'pyg_sphinx_theme'
html_logo = ('https://raw.githubusercontent.com/cgoliver/proteinshake_sphinx_theme/'
             'master/pyg_sphinx_theme/static/img/torch-pdb.png')
html_favicon= ('https://raw.githubusercontent.com/cgoliver/proteinshake_sphinx_theme/'
             'master/pyg_sphinx_theme/static/img/favicon.png')

add_module_names = False
autodoc_member_order = 'bysource'

intersphinx_mapping = {
    'python': ('https://docs.python.org/', None),
    'numpy': ('http://docs.scipy.org/doc/numpy', None),
    'pandas': ('http://pandas.pydata.org/pandas-docs/dev', None),
    'torch': ('https://pytorch.org/docs/master', None),
}


def setup(app):
    def rst_jinja_render(app, _, source):
        rst_context = {'proteinshake': proteinshake}
        source[0] = app.builder.templates.render_string(source[0], rst_context)

    app.connect('source-read', rst_jinja_render)
