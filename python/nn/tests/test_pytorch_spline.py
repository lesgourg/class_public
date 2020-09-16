import pytest
import torch

from ..pytorch_spline import CubicSpline, LinearSpline, digitize


@pytest.fixture
def two_point_spline():
    x = torch.Tensor([0, 1])
    y = torch.Tensor([2, 5])
    return LinearSpline(
            torch.Tensor([0, 1]),
            torch.Tensor([2, 5])
            )


def test_indices():
    xp = torch.Tensor([1, 3, 6])

    x_test = torch.Tensor([0.5, 1, 2, 2.9, 3.1, 4.8, 6, 7])
    indices_expectation = torch.Tensor([-1, 0, 0, 0, 1, 1, 2, 2]).int()

    indices = digitize(xp, x_test)

    assert torch.all(indices == indices_expectation)


def test_single(two_point_spline):
    s = two_point_spline
    assert torch.allclose(s(torch.Tensor([0, 0.5])), torch.Tensor([2, 3.5]))


def test_single_2d_x(two_point_spline):
    s = two_point_spline
    x = torch.Tensor([
        [0, 0.5],
        [0.4, 0.8]
        ])
    y = torch.Tensor([
        [2, 3.5],
        [3.2, 4.4]
        ])
    assert torch.allclose(s(x), y)


def test_single_2d_yp():
    spline = LinearSpline(
            torch.Tensor([0, 1]),
            torch.Tensor([
                [2, 5],
                [4, 1],
            ]))

    x = torch.Tensor([0, 0.5])
    y = torch.Tensor([[2, 3.5], [4, 2.5]])
    y_f = spline(x)
    assert torch.allclose(spline(x), y)


def test_single_2d_x_yp_raises():
    spline = LinearSpline(
            torch.Tensor([0, 1]),
            torch.Tensor([
                [2, 5],
                [4, 1],
            ]))
    x = torch.Tensor([
        [0, 0.5],
        [0.4, 0.8]
        ])

    with pytest.raises(ValueError):
        spline(x)


def test_extrapolate_raises():
    spline = LinearSpline(
            torch.Tensor([1, 3, 6]),
            torch.Tensor([2, 5, 4]),
            )

    with pytest.raises(ValueError):
        spline(torch.Tensor([0]))

    with pytest.raises(ValueError):
        spline(torch.Tensor([7]))



def test_cubic_spline_matches_scipy():
    import numpy as np
    from scipy.interpolate import CubicSpline as ScipyCubicSpline

    x = torch.linspace(-5, 5, 20)
    y = x * torch.cos(x + np.pi / 4) * torch.exp(-x / 3)

    xx = torch.linspace(-5, 5, 100)

    spline = CubicSpline(x, y)
    yy = spline(xx)

    sci_spline = ScipyCubicSpline(x, y, bc_type="natural")
    sci_yy = sci_spline(xx)

    assert torch.allclose(yy, torch.from_numpy(sci_yy).float())
