import torch
import numpy as np
import time


def digitize(boundaries, tensor):
    result = torch.zeros_like(tensor, dtype=torch.int64)
    for boundary in boundaries:
        result += (tensor >= boundary).int()
    return result - 1

class Spline:

    def __init__(self, xp, yp):
        """
        :param xp: x values at which the function has been sampled
        :type xp: torch.Tensor
        :param yp: y values of function evaluated at `x`
        :type yp: torch.Tensor
        """

        self.xp = xp
        self.yp = yp

    def __call__(self, x):
        """Evaluate the spline at the given values `x`.

        :param x: points at which to evaluate the spline
        :type x: torch.Tensor
        :return: values of the spline at `x`.
        :rtype: torch.Tensor
        """

        if not torch.all(x >= self.xp[0]) or not torch.all(x < self.xp[-1]):
            raise ValueError(
                    f"x contains values not in ({self.xp[0]}, {self.xp[-1]})"
                    )

        if not (
                x.ndim == 1 and self.yp.ndim == 1 or
                x.ndim == 2 and self.yp.ndim == 1 or
                x.ndim == 1 and self.yp.ndim == 2):
            raise ValueError("At most one of (self.yp, x) can be 2D!")

        if x.ndim == 1 and self.yp.ndim == 1:
            result_shape = x.shape
        elif x.ndim == 2 and self.yp.ndim == 1:
            result_shape = x.shape
        elif x.ndim == 1 and self.yp.ndim == 2:
            result_shape = self.yp.shape
        else:
            assert False

        result = torch.zeros(result_shape)

        left_boundaries = digitize(self.xp, x)

        return self.interpolate(left_boundaries, x)


class LinearSpline(Spline):

    def interpolate(self, left_boundaries, x):
        right_boundaries = left_boundaries + 1

        xp_lo = self.xp[left_boundaries]
        xp_hi = self.xp[right_boundaries]

        yp_lo = self.yp[..., left_boundaries]
        yp_hi = self.yp[..., right_boundaries]

        result = yp_lo + (yp_hi - yp_lo) / (xp_hi - xp_lo) * (x - xp_lo)

        return result


class DerivativeCubicSpline(Spline):

    def __init__(self, xp, yp, yp_2):
        """
        yp_2: second derivatives
        """

        super().__init__(xp, yp)
        self.yp_2 = yp_2

    def interpolate(self, left_boundaries, x):
        right_boundaries = left_boundaries + 1

        xp_lo = self.xp[left_boundaries]
        xp_hi = self.xp[right_boundaries]

        yp_lo = self.yp[..., left_boundaries]
        yp_hi = self.yp[..., right_boundaries]

        # 2nd derivatives
        m_lo = self.yp_2[..., left_boundaries]
        m_hi = self.yp_2[..., right_boundaries]

        h = xp_hi - xp_lo

        A = (yp_hi - yp_lo) / h - h / 6 * (m_hi - m_lo)
        B = yp_hi - m_hi * h**2 / 6 - A * xp_hi

        return m_lo * (xp_hi - x)**3 / (6 * h) \
                + m_hi * (x - xp_lo)**3 / (6 * h) \
                + A * x \
                + B
                # + (yp_lo - m_lo * h**2 / 6) * (xp_hi - x) / h \
                # + (yp_hi + m_hi * h**2 / 6) * (x - xp_hi) / h
                # + (yp_hi + m_hi * h**2 / 6) * ((x - xp_hi) / h + 1)


def searchsorted(boundaries, tensor):
    scalar = not isinstance(tensor, torch.Tensor)
    if scalar:
        tensor = torch.Tensor([tensor])
    result = torch.zeros_like(tensor, dtype=torch.int64)
    for boundary in boundaries:
        result += (tensor >= boundary).int()
    if scalar:
        return result[0]
    else:
        return result

def diff(x, axis=0, prepend=None):
    high = [slice(None) for _ in range(x.ndim)]
    high[axis] = slice(1, None, None)
    low = [slice(None) for _ in range(x.ndim)]
    low[axis] = slice(None, -1, None)

    result = x[high] - x[low]

    if prepend is not None:
        result = torch.cat((torch.Tensor([prepend]).type(result.dtype), result))

    return result

class CubicSpline:
    def __init__(self, x, y):
        assert x.ndim == 1
        assert y.ndim in (1, 2)

        self.ndim = y.ndim

        if y.ndim == 1:
            y = y[None, :]

        self.x = x.float()
        self.y = y.float()

        self.n = len(self.x) - 1

        self.h = h = diff(self.x)
        self.mu = h[:-1] / (h[1:] + h[:-1])
        self.lamda = 1 - self.mu

        self.build()

    def build(self):
        Q = self.build_matrix()
        v = self.build_differences()

        M = torch.linalg.solve(Q, v)
        self.M = M.to(self.x.device)

    def build_matrix(self):
        Q = torch.zeros((self.n + 1, self.n + 1))

        Q[0, 0] = 1
        Q[-1, -1] = 1

        # Q[0, 0] = 2
        # Q[0, 1] = 1
        # Q[-1, -1] = 2
        # Q[-1, -2] = 1

        for i in range(1, self.n):
            Q[i, i - 1] = self.mu[i - 1]
            Q[i, i] = 2
            Q[i, i + 1] = self.lamda[i - 1]
        return Q

    def build_differences(self):
        result = torch.zeros((self.n + 1, len(self.y)))
        dy = diff(self.y, axis=1).T
        dx = diff(self.x)[:, None]
        rhs = 6 * (dy[1:, :] / dx[1:] - dy[:-1, :] / dx[:-1]) / (self.x[2:, None] - self.x[:-2, None])
        result[1:-1, :] = rhs

        return result

    def __call__(self, x):
        i = searchsorted(self.x, x)
        i[i == self.n + 1] -= 1

        xi = self.x[i]
        xim = self.x[i - 1]
        hi = self.h[i - 1]

        Mi = self.M[i].T
        Mim = self.M[i - 1].T

        dxi = xi - x
        dxim = x - xim

        inv_hi = 1 / hi
        inv_6 = 1 / 6
        inv_6hi = inv_hi * inv_6

        a = dxi**3 * inv_6hi * Mim
        b = dxim**3 * inv_6hi * Mi
        c = (self.y[:, i - 1] * inv_hi - Mim * hi * inv_6) * dxi
        d = (self.y[:, i] * inv_hi - Mi * hi * inv_6) * dxim

        ret = a + b + c + d

        if self.ndim == 1:
            ret = ret.reshape(-1)

        return ret.to(x.device)
