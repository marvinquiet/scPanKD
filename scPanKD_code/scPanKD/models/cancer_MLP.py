import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Function

class GradReverse(Function):
    @staticmethod
    def forward(ctx, x, lambda_):
        ctx.lambda_ = lambda_
        return x.view_as(x)
    
    @staticmethod
    def backward(ctx, grad_output):
        return grad_output.neg() * ctx.lambda_, None

def grad_reverse(x, lambda_):
    return GradReverse.apply(x, lambda_)

class cancer_MLP(nn.Module):
    def __init__(self, dims, input_dim, n_batch, n_classes, dropout=0.5, lambda_adv=0.1):
        super(cancer_MLP, self).__init__()
        self.dims = dims
        self.input_dim = input_dim
        self.n_batch = n_batch
        self.n_classes = n_classes

        layers = [input_dim] + dims
        cancer_MLP_net = []
        for i in range(1, len(layers)):
            if i > 1:
                cancer_MLP_net.append(nn.Dropout(p=dropout))
            linear_layer = nn.Linear(layers[i - 1], layers[i])
            nn.init.xavier_uniform_(linear_layer.weight)
            cancer_MLP_net.append(linear_layer)
            cancer_MLP_net.append(nn.LayerNorm(layers[i]))
            cancer_MLP_net.append(nn.ReLU())
        self.cancer_MLP_net = nn.Sequential(*cancer_MLP_net)

        self.celltype_classifier = nn.Linear(layers[-1], n_classes)
        nn.init.xavier_uniform_(self.celltype_classifier.weight)

        self.cancer_classifier = nn.Linear(layers[-1], n_batch)
        nn.init.xavier_uniform_(self.cancer_classifier.weight)
        # self.batch_predictor = nn.Sequential(
        #     nn.Linear(layers[-1], layers[-1] // 2),
        #     nn.ReLU(),
        #     nn.Linear(layers[-1] // 2, n_batch)
        # )
        self.lambda_adv = lambda_adv

    def forward(self, x, reverse=False):
        features = self.cancer_MLP_net(x)
        if reverse: # domain-invariant feature extraction
            rev_features = grad_reverse(features, self.lambda_adv)
            cancer_pred = self.cancer_classifier(rev_features)
            return cancer_pred
        else:
            class_pred = self.celltype_classifier(features)
            return class_pred, features

    def predict(self, x_test):
        self.eval()
        with torch.no_grad():
            features = self.cancer_MLP_net(x_test)
            logits = self.celltype_classifier(features)
            celltype_prediction = torch.argmax(logits, dim=1)
        return celltype_prediction, features