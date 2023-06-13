
function leaderboard(task) {
    hideDropdown();
    document.getElementById('metric_name').innerHTML = '* Metric: '+metric_names[task];
    document.getElementById('dropdownLabel').innerHTML = task.split('_').map(capitalize).join(' ');
    var board = document.getElementById('leaderboard_table');
    fetch('https://raw.githubusercontent.com/BorgwardtLab/proteinshake/main/leaderboard/'+task+'.json')
        .then((response) => response.json())
        .then((json) => {
            board.innerHTML = '';
            const columns = Object.keys(json[0]);
            var template = ' auto'*columns.length;
            board.style.gridTemplateColumns = new Array(columns.length).fill('auto').join(' ');
            var thead = document.createElement('thead');
            board.appendChild(thead);
            var tr = document.createElement('tr');
            thead.appendChild(tr);
            columns.forEach(x => {
                var th = document.createElement('th');
                th.innerHTML = ['Random Split','Sequence Split','Structure Split'].includes(x) ? x+'*' : x;
                tr.appendChild(th);
            });
            var tbody = document.createElement('tbody');
            board.appendChild(tbody);
            json.sort((a, b) => b['Structure Split'] - a['Structure Split']);
            json.forEach(row => {
                var tr = document.createElement('tr');
                tbody.appendChild(tr);
                columns.forEach(x => {
                    var td = document.createElement('td');
                    td.innerHTML = ['Paper','Code'].includes(x) ? '<a href="'+row[x]+'">Link</a>' : row[x];
                    tr.appendChild(td);
                });
            });
        });
}

function showDropdown() {
    document.getElementById('dropdownContent').style.display = 'block';
}
function hideDropdown() {
    document.getElementById('dropdownContent').style.display = 'none';
}
const capitalize = x => x.charAt(0).toUpperCase() + x.slice(1);

metric_names = {
    'gene_ontology': 'Fmax',
    'enzyme_class': 'Accuracy',
    'protein_family': 'Accuracy',
    'binding_site_detection': 'Matthew\'s correlation coefficient',
    'ligand_affinity': 'Spearman R',
    'protein_protein_interface': 'AUROC (median)',
    'structural_class': 'Accuracy',
    'structure_similarity': 'Spearman R',
};

window.addEventListener('load', ()=>leaderboard('gene_ontology'));
