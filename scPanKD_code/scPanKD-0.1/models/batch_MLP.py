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

class batch_MLP(nn.Module):
    def __init__(self, dims, input_dim, n_batch, n_classes, dropout=0.5, embedding_dim=16, lambda_adv=0.1):
        super(batch_MLP, self).__init__()
        self.dims = dims
        self.input_dim = input_dim
        self.n_batch = n_batch
        self.n_classes = n_classes

        self.batch_embedding = nn.Embedding(n_batch, embedding_dim)
        layers = [input_dim+embedding_dim] + dims
        batch_MLP_net = []
        for i in range(1, len(layers)):
            if i > 1:
                batch_MLP_net.append(nn.Dropout(p=dropout))
            linear_layer = nn.Linear(layers[i - 1], layers[i])
            nn.init.xavier_uniform_(linear_layer.weight)
            batch_MLP_net.append(linear_layer)
            batch_MLP_net.append(nn.LayerNorm(layers[i]))
            batch_MLP_net.append(nn.ReLU())
        self.batch_MLP_net = nn.Sequential(*batch_MLP_net)
        
        self.classifier = nn.Linear(layers[-1], n_classes)
        nn.init.xavier_uniform_(self.classifier.weight)
        
        self.batch_predictor = nn.Sequential(
            nn.Linear(layers[-1], layers[-1] // 2),
            nn.ReLU(),
            nn.Linear(layers[-1] // 2, n_batch)
        )
        self.lambda_adv = lambda_adv

    def forward(self, x, batch_x, reverse=False):
        batch_embed = self.batch_embedding(batch_x)
        features = self.batch_MLP_net(torch.cat((x, batch_embed), dim=1))
        if reverse:
            rev_features = grad_reverse(features, self.lambda_adv)
            batch_pred = self.batch_predictor(rev_features)
            return batch_pred
        else:
            class_pred = self.classifier(features)
            return class_pred, features

    def predict(self, x_test):
        self.eval()
        with torch.no_grad():
            avg_embed = self.batch_embedding.weight.mean(dim=0, keepdim=True)
            batch_vec = avg_embed.expand(x_test.size(0), -1)
            x_input = torch.cat([x_test, batch_vec], dim=1)
            features = self.batch_MLP_net(x_input)
            logits = self.classifier(features)
            probabilities = F.softmax(logits, dim=1)
        return probabilities, features
