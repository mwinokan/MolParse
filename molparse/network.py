class Network(object):
    """Barebones class for visualising a network/graph"""

    def __init__(self, connections):

        import networkx as nx

        self._connections = connections
        self._graph = nx.Graph()

        for u, v in connections:
            self._graph.add_edge(u, v)

    def relax_pos(self, seed=None):
        import networkx as nx

        if seed:
            return nx.spring_layout(self._graph, seed=seed)
        else:
            return nx.spring_layout(self._graph)

    def plot(self, seed=None, show=True):

        import networkx as nx
        import matplotlib.pyplot as plt

        pos = self.relax_pos(seed)

        options = {
            "font_size": 36,
            "node_size": 3000,
            "node_color": "white",
            "edgecolors": "black",
            "linewidths": 5,
            "width": 5,
        }

        # draw network
        nx.draw_networkx(self._graph, pos, **options)

        ax = plt.gca()
        ax.margins(0.2)
        plt.axis("off")

        if show:
            plt.show()

        return ax
