var representation = 'graph';
var framework = 'torch';

const lines = {
    import: {
        pyg: 'from torch_geometric.loader import DataLoader',
        dgl: 'from torch.utils.data import DataLoader',
        nx: '',
        torch: 'from torch.utils.data import DataLoader',
        tf: 'from tf.data import Dataset',
        np: '',
    },
    comment_1: {
        graph: '# Convert them to graphs with an epsilon neighborhood of 8 Angstrom',
        voxel: '# Convert them to voxels with a voxelsize of 10 Angstrom',
        point: '# Convert them to point clouds',
    },
    comment_2: {
        pyg: '# Load into PyTorch-Geometric data structures',
        dgl: '# Load into DGL data structures',
        nx: '# Load into NetworkX data structures',
        torch: '# Load into PyTorch data structures',
        tf: '# Load into Tensorflow data structures',
        np: '# Load into Numpy data structures',
    },
    convert: {
        graph: 'graph(eps=8)',
        voxel: 'voxel(voxelsize=10)',
        point: 'point()',
    },
    comment_3: {
        pyg: '# Training using native data loaders',
        dgl: '# Training using native data loaders',
        nx: '# Training',
        torch: '# Training using native data loaders',
        tf: '# Training using native data loaders',
        np: '# Training',
    },
    loader: {
        pyg: 'DataLoader(task.train)',
        dgl: 'DataLoader(task.train)',
        nx: 'task.train',
        torch: 'DataLoader(task.train)',
        tf: 'Dataset.from_tensor_slices(task.train)',
        np: 'task.train',
    },
};

function quickstart() {
    var representations = document.getElementById('representations');
    [...representations.children].forEach(x=>x.classList.remove('selected'));
    var selected = document.getElementById(representation);
    selected.classList.add('selected');

    if (representation == 'graph') {
        ['torch','tf','np'].forEach(x=>document.getElementById(x).style.display='none');
        ['pyg','dgl','nx'].forEach(x=>document.getElementById(x).style.display='flex');
        if (['torch','tf','np'].includes(framework)) framework = {torch:'pyg',tf:'dgl',np:'nx'}[framework];
    } else {
        ['torch','tf','np'].forEach(x=>document.getElementById(x).style.display='flex');
        ['pyg','dgl','nx'].forEach(x=>document.getElementById(x).style.display='none');
        if (['pyg','dgl','nx'].includes(framework)) framework = {pyg:'torch',dgl:'tf',nx:'np'}[framework];
    }

    var frameworks = document.getElementById('frameworks');
    [...frameworks.children].forEach(x=>x.classList.remove('selected'));
    var selected = document.getElementById(framework);
    selected.classList.add('selected');

    var code = document.getElementById('code');
    code.innerHTML =
`from proteinshake.tasks import GeneOntologyTask
${lines.import[framework]}

# Use proteins with Gene Ontology annotations
${lines.comment_1[representation]}
${lines.comment_2[framework]}
task = GeneOntologyTask().to_${lines.convert[representation]}.${framework}()

${lines.comment_3[framework]}
for batch in ${lines.loader[framework]}:
    model.train_step(batch) # your model training goes here

# Evaluation with the provided metrics
prediction = model.inference(task.test)
metrics = task.evaluate(task.test_targets, prediction)

print(metrics)`;
    Prism.highlightAll();
}


var task = 'gene_ontology';
function leaderboard() {
    hideDropdown();
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
                th.innerHTML = x;
                tr.appendChild(th);
            });
            var tbody = document.createElement('tbody');
            board.appendChild(tbody);
            json.sort((a, b) => a['Structure Split'].localeCompare(b['Structure Split']));
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
function copyToClipboard() {
    navigator.clipboard.writeText('pip install proteinshake');
    var copy = document.getElementById('copy');
    copy.src = 'images/done.svg';
    setTimeout(x=>copy.src = 'images/copy.svg', 1000);
}
const capitalize = x => x.charAt(0).toUpperCase() + x.slice(1);


window.addEventListener('load', quickstart);
window.addEventListener('load', leaderboard);
