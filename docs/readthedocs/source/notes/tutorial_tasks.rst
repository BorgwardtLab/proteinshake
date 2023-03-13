Training a protein structure model
===================================

Now that we understand how to create datasets and cast them to an ML framework we can take a look at the tasks class.
The ``proteinshake.Tasks`` class exposes two main functionalities for various tasks: dataset splitting, and model evaluation.
This, along with the ProteinShake datasets, leaves the user to focus only on coding the model.


Basic Usage
~~~~~~~~~~~~

The user interaction with the task object is quite simple. At a glance it looks like this::

        >>> from proteinshake.tasks import EnzymeCommissionTask
        >>> dataset = task.dataset.to_graph(eps=8).pyg()
        >>> pred = model(dataset[task.train_index]) # assuming you have implemented model() elsewhere
        >>> task.evaluate(pred)
        {'precision': 0.5333515066547034, 'recall': 0.4799021029676011, 'accuracy': 0.6675514266755143}



Next we provide a full tutorial including model construction on the Enzyme Commission task.

Predicting enzyme class with ProteinShake 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: https://colab.research.google.com/assets/colab-badge.svg 
   :target: https://colab.research.google.com/github/BorgwardtLab/proteinshake/blob/main/examples/sup_enzyme_commission_with_gnn.ipynb 

We will use a simple GNN model, namely GCN, and evaluate its performance for enzyme commission prediction.
The model can be trained with either CPU or GPU, but GPU is recommended for faster computation.



Let's get some imports::


        import copy
        from tqdm import tqdm
        import torch
        from torch import nn
        import torch.nn.functional as F
        from proteinshake import tasks as ps_tasks



We import the task object and convert the protein 3D structures to graphs using an epsilon cutoff of 8 Angstroms::

        def transform(data):
            data, protein_dict = data
            data.y = task.target(protein_dict)
            return data


        datapath = './data/ec'
        task = ps_tasks.EnzymeCommissionTask(root=datapath)
        dset = task.dataset
                dset = dset.to_graph(eps=8.0).pyg(
            transform=transform
        )

We can now create data loaders for train/val/test sets provided by ProteinShake::

        from torch.utils.data import Subset
        from torch_geometric.loader import DataLoader
        batch_size = 100
        train_loader = DataLoader(Subset(dset, task.train_index), batch_size=batch_size,
                                  shuffle=True, num_workers=0)
        val_loader = DataLoader(Subset(dset, task.val_index), batch_size=batch_size,
                                shuffle=False, num_workers=0)
        test_loader = DataLoader(Subset(dset, task.test_index), batch_size=batch_size,
                                 shuffle=False, num_workers=0)



Next, we build a simple graph convolution network using pytorch geometric.


This is the convolution layer which is applied to messages from neighoring residues::

        import torch_geometric.nn as gnn
        from torch_geometric import utils
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



Now the graph neural network::

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


We build a GCN model with 5 layers and 64 hidden dimensions::


        embed_dim = 64
        num_layers = 5

        model = GNN_graphpred(
            task.num_classes,
            embed_dim,
            num_layers,
        )


Build an optimizer and define the train and test function::

        lr = 0.001
        optimizer = torch.optim.AdamW(
            model.parameters(),
            lr=lr
        )

        criterion = nn.CrossEntropyLoss()
        # set device
        device = torch.device(torch.cuda.current_device()) \
        if torch.cuda.is_available() else torch.device('cpu')



Define the logic to apply at each training epoch::

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


ProteinShake provides an evaluation function for each task ``task.evaluate(y_pred)``::


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


Now the training stage::

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


Let's see the model performance on the different evaluation metrics provided::

        >>> test_scores = eval_epoch(model, test_loader)
        >>> test_scores
        {'precision': 0.5333515066547034, 'recall': 0.4799021029676011, 'accuracy': 0.6675514266755143}


:download:`Download source code for this example. <code/ec_task.py>`
