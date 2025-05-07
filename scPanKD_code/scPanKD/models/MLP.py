import torch
import torch.nn as nn
import torch.nn.functional as F

class MLP(nn.Module):
    def __init__(self, dims, input_shape, n_classes, random_seed=1993):
        super(MLP, self).__init__()
        torch.manual_seed(random_seed)
        self.layers = nn.ModuleList()
        self.input_shape = input_shape
        self.n_classes = n_classes

        # Define layers
        prev_dim = input_shape
        for dim in dims:
            self.layers.append(nn.Sequential(
                nn.Dropout(p=0.5),
                nn.Linear(prev_dim, dim),
                nn.ReLU()
            ))
            prev_dim = dim
        self.output_layer = nn.Linear(prev_dim, n_classes)

    def forward(self, x):
        for layer in self.layers:
            x = layer(x)
        x = self.output_layer(x)
        return x

    def compile(self, optimizer='adamW', learning_rate=1e-4):
        if optimizer == 'adam':
            self.optimizer = torch.optim.Adam(filter(lambda p: p.requires_grad, self.parameters()), lr=learning_rate)
        if optimizer == 'adamW':
            self.optimizer = torch.optim.AdamW(filter(lambda p: p.requires_grad, self.parameters()), lr=learning_rate)
        self.criterion = nn.CrossEntropyLoss()

    def fit(self, x_train, y_train, batch_size=16, max_epochs=100):
        dataset = torch.utils.data.TensorDataset(x_train, y_train)
        dataloader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=True)
        for epoch in range(max_epochs):
            self.train()
            for batch_x, batch_y in dataloader:
                self.optimizer.zero_grad()
                outputs = self.forward(batch_x)
                loss = self.criterion(outputs, batch_y)
                loss.backward()
                self.optimizer.step()

    def predict(self, x_test):
        self.eval()
        with torch.no_grad():
            outputs = self.forward(x_test)
            probabilities = F.softmax(outputs, dim=1)
        return probabilities


class batch_MLP(nn.Module):
    def __init__(self, dims, input_dim, n_batch, n_classes, dropout=0.5, random_seed=1993):
        super(batch_MLP, self).__init__()
        torch.manual_seed(random_seed)
        self.dims = dims
        self.input_dim = input_dim
        self.n_batch = n_batch
        self.n_classes = n_classes

        layers = []
        batch_layers = []
        # Input layers
        layers.append(nn.Sequential(
            nn.Dropout(p=dropout_rate),
            nn.Linear(self.input_shape, self.dims[0]),
            nn.ReLU()
        ))
        batch_layers.append(nn.Linear(self.n_batch, self.dims[0]))
        # Hidden layers
        for i in range(len(self.dims) - 1):
            layers.append(nn.Sequential(
                nn.Dropout(p=dropout_rate),
                nn.Linear(self.dims[i], self.dims[i + 1]),
                nn.ReLU()
            ))
            batch_layers.append(nn.Linear(self.n_batch, self.dims[i + 1]))

        # Output layer
        layers.append(nn.Linear(self.dims[-1], self.n_classes))
        self.layers = nn.ModuleList(layers)
        self.batch_layers = nn.ModuleList(batch_layers)

    def forward(self, x, batch_x):
        for i in range(len(self.layers) - 1):
            x = self.layers[i](x)
            batch_x_proj = self.batch_layers[i](batch_x)
            x = x + batch_x_proj
        x = self.layers[-1](x)
        return x

    def compile(self, optimizer='adamW', learning_rate=1e-4):
        if optimizer == 'adam':
            self.optimizer = torch.optim.Adam(filter(lambda p: p.requires_grad, self.parameters()), lr=learning_rate)
        if optimizer == 'adamW':
            self.optimizer = torch.optim.AdamW(filter(lambda p: p.requires_grad, self.parameters()), lr=learning_rate)
        self.criterion = nn.CrossEntropyLoss()

    def fit(self, x_train, b_train, y_train, batch_size=16, max_epochs=100):
        dataset = torch.utils.data.TensorDataset(x_train, b_train, y_train)
        dataloader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=True)

        for epoch in range(max_epochs):
            self.train()
            for batch_x, batch_b, batch_y in dataloader:
                self.optimizer.zero_grad()
                outputs = self.forward(batch_x, batch_b)
                loss = self.criterion(outputs, batch_y)
                loss.backward()
                self.optimizer.step()

    def predict(self, x_test, b_test):
        self.eval()
        with torch.no_grad():
            outputs = self.forward(x_test, b_test)
            probabilities = F.softmax(outputs, dim=1)
        return probabilities
