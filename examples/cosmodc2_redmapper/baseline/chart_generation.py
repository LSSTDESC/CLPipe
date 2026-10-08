import sys
import os

sys.path.insert(0, "/sps/lsst/groups/clusters/cl_pipeline_project/TXPipe")
sys.path.insert(0, "/sps/lsst/ebarroso/CLPipe")
from ceci import Pipeline

config_path = "full_pipeline_chart.yml"
pipe_config = Pipeline.build_config(config_path, dry_run=True)

pipeline = Pipeline.create(pipe_config)

graph = pipeline.make_flow_chart(filename=None)

graph.graph_attr.update(dpi="300", fontsize="30", ranksep="0.6", nodesep="0.35")
graph.node_attr.update(fontsize="18", fontname="Helvetica", margin="0.15,0.08")
graph.edge_attr.update(fontsize="12", fontname="Helvetica")

for node in graph.nodes_iter():
    node.attr["fontsize"] = "24"

out_path = "flowchart_pipeline.pdf"
graph.draw(out_path, prog="dot")
print(f"Saved {out_path}")
