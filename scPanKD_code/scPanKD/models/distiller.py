## from Keras Knowledge Distillation: https://keras.io/examples/vision/knowledge_distillation/
## for predicting scATAC-seq using another scATAC-seq
import torch
import torch.nn as nn
import torch.nn.functional as F

class Distiller(nn.Module):
    def __init__(self, student, teacher):
        super(Distiller, self).__init__()
        self.teacher = teacher
        self.student = student

    def compile(self, optimizer, student_loss_fn, distillation_loss_fn, alpha=0.1, temperature=3):
        """
        Configure the distiller.

        Args:
            optimizer: PyTorch optimizer for the student weights
            student_loss_fn: Loss function of difference between student
                predictions and ground-truth
            distillation_loss_fn: Loss function of difference between soft
                student predictions and soft teacher predictions
            alpha: weight to student_loss_fn and 1-alpha to distillation_loss_fn
            temperature: Temperature for softening probability distributions.
                Larger temperature gives softer distributions.
        """
        self.optimizer = optimizer
        self.student_loss_fn = student_loss_fn
        self.distillation_loss_fn = distillation_loss_fn
        self.alpha = alpha
        self.temperature = temperature

    def train_step(self, x, y):
        # Forward pass of teacher
        with torch.no_grad():
            teacher_predictions = self.teacher(x)

        # Forward pass of student
        student_predictions = self.student(x)

        # Compute losses
        student_loss = self.student_loss_fn(student_predictions, y)
        distillation_loss = self.distillation_loss_fn(
            F.softmax(teacher_predictions / self.temperature, dim=1),
            F.softmax(student_predictions / self.temperature, dim=1),
        )
        loss = self.alpha * student_loss + (1 - self.alpha) * distillation_loss

        # Backward pass and optimization
        self.optimizer.zero_grad()
        loss.backward()
        self.optimizer.step()

        return {
            "student_loss": student_loss.item(),
            "distillation_loss": distillation_loss.item(),
            "total_loss": loss.item()
        }

    def test_step(self, x, y):
        # Compute predictions
        with torch.no_grad():
            student_predictions = self.student(x)

        # Calculate the loss
        student_loss = self.student_loss_fn(student_predictions, y)

        return {
            "student_loss": student_loss.item()
        }

