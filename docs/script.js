var representation = 'graph';
var framework = 'torch';

/*
from proteinshake.tasks import GeneOntologyTask
from torch_geometric.data import DataLoader

# Use proteins with Gene Ontology annotations
# Convert them to graphs with an epsilon neighborhood of 8 Angstrom
# Load into PyTorch-Geometric data structures
task = GeneOntologyTask().to_graph(eps=8).pyg()

# Training using native data loaders
for X,y in DataLoader(task.X_train, task.y_train, batchsize=32):
    model.train_step(X,y) # your model training goes here

# Evaluation with the provided metrics
prediction = model.inference(task.X_test)
metrics = task.eval(prediction)

print(metrics)
*/

function quickstart() {
    var representations = document.getElementById('representations');
    [...representations.children].forEach(x=>x.classList.remove('selected'));
    var selected = document.getElementById(representation);
    selected.classList.add('selected');

    var frameworks = document.getElementById('frameworks');
    [...frameworks.children].forEach(x=>x.classList.remove('selected'));
    var selected = document.getElementById(framework);
    selected.classList.add('selected');

    var code = document.getElementById('code');
    code.innerHTML = `
    from proteinshake.tasks import GeneOntologyTask
    from torch_geometric.data import DataLoader

    # Use proteins with Gene Ontology annotations
    # Convert them to graphs with an epsilon neighborhood of 8 Angstrom
    # Load into PyTorch-Geometric data structures
    task = GeneOntologyTask().to_graph(eps=8).pyg()

    # Training using native data loaders
    for X,y in DataLoader(task.X_train, task.y_train, batchsize=32):
        model.train_step(X,y) # your model training goes here

    # Evaluation with the provided metrics
    prediction = model.inference(task.X_test)
    metrics = task.eval(prediction)

    print(metrics)
    `;
    Prism.highlightAll();
}


window.addEventListener('load', quickstart);
