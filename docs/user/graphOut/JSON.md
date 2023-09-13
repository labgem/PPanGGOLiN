The json's file content corresponds to the .gexf file content, but in json rather than gexf file format. It follows the 'node-link' format as shown in [this example](https://observablehq.com/@d3/force-directed-graph) in javascript, or as used in the [networkx](https://networkx.github.io/documentation/stable/reference/readwrite/json_graph.html) python library and it should be usable with both [D3js](https://d3js.org/) and [networkx](https://networkx.github.io/documentation/stable/index.html), or any other software or library that supports this format.

The json can be generated using the 'write' subcommand as such : 

`ppanggolin write -p pangenome.h5 --json`