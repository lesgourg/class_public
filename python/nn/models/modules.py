import torch
from torch import nn


class DenseNet(nn.Module):

    def __init__(self, n_inputs, n_outputs, hidden_units):
        super().__init__()

        self.net = self._build_net(n_inputs, n_outputs, hidden_units)

    def _build_net(self, n_in, n_out, hidden):
        if isinstance(hidden, int):
            hidden = [hidden]

        layers = []
        last = n_in
        for h in hidden:
            layers.append(nn.Linear(last, h))
            layers.append(nn.ReLU(True))
            last = h
        layers.append(nn.Linear(last, n_out))
        return nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)


class DampingNet(nn.Module):

    def __init__(self, n_inputs, n_outputs, hidden_units):
        super().__init__()

        self.dense = DenseNet(n_inputs, n_outputs, hidden_units)

    def forward(self, x, k, k_d):
        factors = self.dense(x)
        envelopes = torch.exp(-torch.pow(
            k[None, None, :] / k_d[None, :, None] * (1 + factors.T[..., None]),
            2))
        return envelopes


class PhaseNet(nn.Module):

    def __init__(self, n_inputs, n_phases, max_order, hidden_units, zero_phases=None):
        """
        zero_phases: [(phase_index, coefficient_index), ...]
        """
        super().__init__()

        self.n_phases = n_phases
        self.max_order = max_order
        n_outputs = n_phases * (max_order + 1)

        self.zero_phases = zero_phases if zero_phases is not None else []

        self.dense = DenseNet(n_inputs, n_outputs, hidden_units)

    def forward(self, x, k, r_s):
        # shape: (n_phases, n_coefficients, n_tau)
        values = self.dense(x).reshape(self.n_phases, self.max_order + 1, -1)

        # (n_tau, n_k)
        arg = k[None, :] * r_s[:, None]

        # (n_phases, n_coefficients, n_tau, n_k)
        powers = torch.arange(self.max_order + 1)[None, :, None, None].to(x.device)

        # (n_phases, n_coefficients, n_tau, n_k)
        multiplied = values[..., None] * torch.pow(arg, powers)

        for i_phase, i_c in self.zero_phases:
            multiplied[i_phase][i_c] = 0

        # add usual k dep. phase
        if self.max_order >= 1:
            multiplied[1] += arg

        # (n_phases, n_tau, n_k)
        phases = multiplied.sum(dim=1)

        return phases


