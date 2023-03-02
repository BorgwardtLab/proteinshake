import copy
from tqdm import tqdm
import torch
from torch import nn
from torch.utils.data import Subset
from torch_geometric.loader import DataLoader
import torch_geometric.nn as gnn
from torch_geometric import utils
import torch.nn.functional as F
from proteinshake import tasks as ps_tasks


def transform(data):
    data, protein_dict = data
    data.y = task.target(protein_dict)
    return data

class GCNConv(gnn.MessagePassing):
    def __init__(self, embed_dim=256, use_edge_attr=False):
        super().__init__(aggr='add')
        self.use_edge_attr = use_edge_attr

        self.linear = nn.Linear(embed_dim, embed_dim)
        self.root_emb = nn.Embedding(1, embed_dim)
        self.edge_encoder = nn.Linear(embed_dim, embed_dim)

    def forward(self, x, edge_index, edge_attr=None):
        x = self.linear(x)
        if self.use_edge_attr and edge_attr is not None:
            edge_attr = self.edge_encoder(edge_attr)

        row, col = edge_index

        deg = utils.degree(row, x.size(0), dtype = x.dtype) + 1
        deg_inv_sqrt = deg.pow(-0.5)
        deg_inv_sqrt[deg_inv_sqrt == float('inf')] = 0

        norm = deg_inv_sqrt[row] * deg_inv_sqrt[col]

        return self.propagate(
            edge_index, x=x, edge_attr = edge_attr, norm=norm) + F.relu(
            x + self.root_emb.weight) * 1./deg.view(-1,1)

    def message(self, x_j, edge_attr, norm):
        return norm.view(-1, 1) * F.relu(x_j + edge_attr)


class GNN(nn.Module):
        def __init__(self, embed_dim=256, num_layers=3, dropout=0.0,
                      use_edge_attr=False):
                super().__init__()
                self.embed_dim = embed_dim
                self.num_layers = num_layers
                self.dropout = dropout

                self.x_embedding = nn.Embedding(20, embed_dim)

                gnn_model = GCNConv
                self.gnns = nn.ModuleList()
                for _ in range(num_layers):
                        self.gnns.append(gnn_model(embed_dim, use_edge_attr=use_edge_attr))

                self.batch_norms = nn.ModuleList()
                for _ in range(num_layers):
                        self.batch_norms.append(nn.BatchNorm1d(embed_dim))

        def forward(self, data):
                x, edge_index, edge_attr = data.x, data.edge_index, data.edge_attr

                output = self.x_embedding(x)

                for layer in range(self.num_layers):
                    output = self.gnns[layer](output, edge_index, edge_attr)
                    output = self.batch_norms[layer](output)

                    if layer == self.num_layers - 1:
                        output = F.dropout(output, self.dropout, training=self.training)
                    else:
                        output = F.dropout(F.relu(output), self.dropout, training=self.training)

                return output

@torch.no_grad()
def eval_epoch(model, loader):
    model.eval()

    y_true = []
    y_pred = []

    for step, batch in enumerate(loader):
        batch = batch.to(device)
        y_hat = model(batch)

        y_true.append(batch.y.cpu())
        y_pred.append(y_hat.cpu())

    y_true = torch.cat(y_true, dim = 0).numpy()
    y_pred = torch.vstack(y_pred).numpy()
    y_pred = y_pred.argmax(-1)
    scores = task.evaluate(y_true, y_pred)
    return scores

def train_epoch(model):
    model.train()

    running_loss = 0.
    for step, batch in enumerate(train_loader):
        size = len(batch.y)
        batch = batch.to(device)

        optimizer.zero_grad()
        y_hat = model(batch)

        loss = criterion(y_hat, batch.y)

        loss.backward()
        optimizer.step()

        running_loss += loss.item() * size

    n_sample = len(train_loader.dataset)
    epoch_loss = running_loss / n_sample
    return epoch_loss


if __name__ == "__main__":
    datapath = './data/ec'
    task = ps_tasks.EnzymeClassTask(root=datapath)
    dset = task.dataset

    dset = dset.to_graph(eps=8.0).pyg(
    transform=transform
    )

    batch_size = 100
    embed_dim = 64
    num_layers = 5

    criterion = nn.CrossEntropyLoss()
    lr = 0.001

    # set device
    device = torch.device(torch.cuda.current_device()) \
    if torch.cuda.is_available() else torch.device('cpu')

    train_loader = DataLoader(Subset(dset, task.train_index), batch_size=batch_size,
                              shuffle=True, num_workers=0)
    val_loader = DataLoader(Subset(dset, task.val_index), batch_size=batch_size,
                            shuffle=False, num_workers=0)
    test_loader = DataLoader(Subset(dset, task.test_index), batch_size=batch_size,
                         shuffle=False, num_workers=0)



    model = GNN_graphpred(
        task.num_classes,
        embed_dim,
        num_layers,
    )


    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=lr
    )



    model.to(device)
    epochs = 20 # we train only 20 epochs here, but more epochs may result in better performance.

    best_val_score = 0.0
    pbar = tqdm(range(epochs))
    for epoch in pbar:
        train_loss = train_epoch(model)
        val_scores = eval_epoch(model, val_loader)
        val_score = val_scores['accuracy']
        postfix = {'train_loss': train_loss, 'val_acc': val_score}
        pbar.set_postfix(postfix)

        if val_score > best_val_score:
            best_val_score = val_score
            best_weights = copy.deepcopy(model.state_dict())

    model.load_state_dict(best_weights)


    test_scores = eval_epoch(model, test_loader)
    print(test_scores)
