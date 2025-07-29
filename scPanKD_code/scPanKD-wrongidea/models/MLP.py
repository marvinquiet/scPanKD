import torch
import torch.nn as nn
import torch.nn.functional as F

class MLP(nn.Module):
    def __init__(self, dims, input_dim, n_classes, dropout=0.5):
        super(MLP, self).__init__()
        self.dims = dims
        self.input_dim = input_dim
        self.n_classes = n_classes

        layers = [input_dim] + dims
        MLP_net = []
        for i in range(1, len(layers)):
            MLP_net.append(nn.Dropout(p=dropout))
            linear_layer = nn.Linear(layers[i - 1], layers[i])
            nn.init.xavier_uniform_(linear_layer.weight)
            MLP_net.append(linear_layer)
            MLP_net.append(nn.LayerNorm(layers[i]))
            MLP_net.append(nn.ReLU())
        linear_layer = nn.Linear(layers[-1], n_classes)
        nn.init.xavier_uniform_(linear_layer.weight)
        MLP_net.append(linear_layer)
        self.MLP_net = nn.Sequential(*MLP_net)

    def forward(self, x):
        return self.MLP_net(x)

    def predict(self, x_test):
        self.eval()
        with torch.no_grad():
            outputs = self.forward(x_test)
            probabilities = F.softmax(outputs, dim=1)
        return probabilities