



Adjacency = [
0 1 1 0 1 0 0 1 0 0 0 1 0 1 0 0 0 0;
1 0 0 0 0 0 1 0 0 0 0 1 0 1 1 0 0 0;
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
0 0 0 0 1 1 0 0 0 1 1 0 0 0 0 0 0 0;
1 0 0 1 0 0 0 0 1 0 1 0 1 0 0 0 0 0;
0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0;
1 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0;
0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0;
0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
0 0 0 1 1 0 0 0 0 0 0 0 0 1 0 0 0 0;
1 1 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0;
0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1;
1 1 0 0 0 0 0 0 0 0 1 1 0 0 1 0 0 0;
0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0;
]
#n_verts = 18

edges  = edge_list(Adjacency)

colours = [1, 1, 2, 3, 3, 3, 1, 3, 3, 3, 3, 1, 2, 1, 1, 1, 2, 2]

Model = β_mod(Adjacency, colours)
#The number of variables is 18 + 3 for the vertices and all pairs of distinct colours.
#There are 28 functions, the number of edges.
#The toric parametrisation is finite since rank(Model) == 21.

θ, q = to_system(Model)
