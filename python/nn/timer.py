import time
import contextlib

class Timer:
    """
    Simple & convenient timing class.
    """

    def __init__(self):
        self._times = {}
        self._start = {}

    @property
    def times(self):
        return self._times

    def start(self, name):
        self._start[name] = time.perf_counter()

    def end(self, name):
        if not name in self._start:
            raise ValueError(f"Cannot end time '{name}' because it is not started!")
        self._times[name] = time.perf_counter() - self._start[name]
        del self._start[name]

    @contextlib.contextmanager
    def __call__(self, name):
        self.start(name)
        yield
        self.end(name)

