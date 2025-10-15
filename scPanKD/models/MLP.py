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
            linear_layer = nn.Linear(layers[i - 1], layers[i])
            nn.init.xavier_uniform_(linear_layer.weight)
            MLP_net.append(linear_layer)
            MLP_net.append(nn.LayerNorm(layers[i]))
            MLP_net.append(nn.ReLU())
        self.MLP_net = nn.Sequential(*MLP_net)
        self.celltype_classifier = nn.Linear(layers[-1], n_classes)
        nn.init.xavier_uniform_(self.celltype_classifier.weight)

    def forward(self, x):
        features = self.MLP_net(x)
        class_pred = self.celltype_classifier(features)
        return class_pred

    def predict(self, x_test):
        self.eval()
        with torch.no_grad():
            features = self.MLP_net(x_test)
            logits = self.celltype_classifier(features)
            probabilities = F.softmax(logits, dim=1)
        return probabilities, features