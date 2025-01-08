from ppanggolin.genome import Feature
from ppanggolin.edge import Edge

class Intergenomic(Feature):
    """Constructor method
        :param gene_id: Identifier of the gene
    """
    def __init__(self, edge: Edge): # check if should add ID to Edge
        self.edge = edge
        self.sequence = None
        self.source = self.edge.source
        self.target = self.edge.target

    @property

    def add_edge(self, edge: Edge):
        if not isinstance(edge, str):
            raise TypeError(
                f"'str' type was expected but you provided a '{type(edge)}' type object"
            )
        self.edge = edge

    @prperty
    def extract_sequence(self):
        if self.sequence is None:







