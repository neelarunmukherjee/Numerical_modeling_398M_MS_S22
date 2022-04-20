mu = 1; % viscosity - value does not matter here

% Definition of pressure grid
Gridp.xmin = 0; Gridp.xmax = 1; Gridp.Nx = 5;
Gridp.ymin = 0; Gridp.ymax = 1; Gridp.Ny = 5;

Grid = build_stokes_grid(Gridp);