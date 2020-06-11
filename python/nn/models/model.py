"""
Abstract Base Class for all network models.
"""

from abc import ABC, abstractmethod

import torch
import torch.nn as nn

class Model(ABC, nn.Module):
    def __init__(self, k):
        self.k = k

    def criterion(self):
        return nn.MSELoss()

    def optimizer(self):
        return torch.optim.Adam(self.parameters(), lr=self.learning_rate)

    def tau_training(self):
        """
        Returns the tau array on which to train the model.
        This method may return None if the model should be trained
        on the tau array obtained from CLASS during generation of training
        data (and does so by default).
        """
        return None

    @abstractmethod
    def lr_scheduler(self, optimizer):
        """
        Given an `optimizer`, returns a LR scheduler.
        """
        pass

    @abstractmethod
    def epochs(self):
        """Returns the number of epochs to train the model for."""
        pass

    @abstractmethod
    def required_inputs(self):
        """
        Returns a list of input names to be passed to this
        model during training.
        """
        pass

    @abstractmethod
    def source_functions(self):
        """
        Returns a list of source function names to be passed
        to this model during training.
        """
        pass
