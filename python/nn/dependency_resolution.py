import functools
import inspect

class Graph:

    def __init__(self):
        self.graph = {}

    def require(self, *quantities):
        quantity_names = [q if isinstance(q, str) else q.__name__ for q in quantities]
        def decorator(func):
            self.graph[func.__name__] = (quantity_names, func)
            return func
        return decorator

    def register(self, func):
        signature = inspect.signature(func)
        quantities = list(signature.parameters)
        self.require(*quantities)(func)

    def resolve(self, node, stack=None, visited=None):
        if stack is None:
            stack = []
        if visited is None:
            visited = {node: False for node in self.graph}

        visited[node] = True

        adjacent = self.graph[node][0]
        for adj in adjacent:
            if adj in self.graph and not visited[adj]:
                self.resolve(adj, stack, visited)

        stack.append(node)
        return stack

    def resolve_all(self):
        return {node: self.resolve(node) for node in self.graph}

    def evaluator(self):
        return Evaluator(self, self.resolve_all())

class Evaluator:

    def __init__(self, graph, resolution_chains):
        self.graph = graph
        self.resolution_chains = resolution_chains

    def evaluate(self, qty, cache=None):
        if cache is None:
            cache = {}

        # If requested quantity has already been computed earlier OR 
        # was in the cache initially, simply return it
        if qty in cache:
            return cache[qty]

        # Otherwise, walk the dependencies of requested quantity
        # in topological order as determind by graph 
        for dep in self.resolution_chains[qty]:
            # If the current dependency has not been calculated yet (i.e. is not 
            # in the cache), compute it and store it to the cache
            if dep not in cache:
                # Obtain all the arguments for the current quantity.
                # These are guaranteed to already be in cache, because
                # self.resolution_chains[...] are topologically sorted,
                # i.e. if B depends on A, A will already be in the cache
                args = [cache[ddep] for ddep in self.graph.graph[dep][0]]
                # print("dep: {}; args: {}".format(dep, args))
                func = self.graph.graph[dep][1]
                cache[dep] = func(*args)

        # The last entry in a resolution_chain is the requested quantity itself,
        # e.g. self.resolution_chains[X][-1] == X.
        # Hence, this quantity will be in cache as well and can simply be returned
        return cache[qty]

    def evaluate_many(self, quantities, cache=None):
        if cache is None:
            cache = {}
        return {quantity: self.evaluate(quantity, cache=cache) for quantity in quantities}

    def evaluate_all(self, cache=None):
        return self.evaluate_many(self.graph.graph.keys(), cache=cache)

