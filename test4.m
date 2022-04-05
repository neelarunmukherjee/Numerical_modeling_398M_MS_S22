Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 4;
Grid.ymin = 0; Grid.ymax = 1; Grid.Ny = 3;
Grid = build_grid2(Grid);
[D,G,C,I,M]=build_ops2(Grid);

subplot 121
spy(D), title 'D'
subplot 122
spy(G), title 'G'