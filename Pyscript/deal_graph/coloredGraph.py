node_colours = []
for i in range(node_count):
    if i in efaecalis_list:
        node_colours.append("red")
    elif i in saureus_list:
        node_colours.append("green")
    else:
        node_colours.append("white")
out_fig_name = "coloured_assembly_graph.png"
g.vs["color"] = node_colours
visual_style = {}
# Set bbox and margin
visual_style["bbox"] = (1500,1500)
visual_style["margin"] = 30
# Set vertex size
visual_style["vertex_size"] = 35
# Set vertex lable size
visual_style["vertex_label_size"] = 15
# Don't curve the edges
visual_style["edge_curved"] = False
# Set the layout
visual_style["layout"] = my_layout
# Plot the graph
plot(g, out_fig_name, **visual_style)