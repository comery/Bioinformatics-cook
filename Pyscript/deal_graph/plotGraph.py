out_fig_name = "assembly_graph.png"
visual_style = {}
# Set bbox and margin
visual_style["bbox"] = (1500,1500)
visual_style["margin"] = 30
# Set vertex colours
visual_style["vertex_color"] = 'white'
# Set vertex size
visual_style["vertex_size"] = 35
# Set vertex lable size
visual_style["vertex_label_size"] = 15
# Don't curve the edges
visual_style["edge_curved"] = False
# Set the layout
my_layout = g.layout_fruchterman_reingold()
visual_style["layout"] = my_layout
# Plot the graph
plot(g, out_fig_name, **visual_style)