function [D,G,C,I,M]=build_ops(Grid)
% author: 
% date: 
% description:
% This function computes the discrete operator matrices on a
% staggered grid using central difference approximations. The
% discrete gradient assumes homogeneous boundary conditions.
% Input:
% Grid = structure containing all pertinent information about the grid.
% Output:
% D = Nx by Nx+1 discrete divergence matrix 
% G = Nx+1 by Nx discrete gradient matrix
% C = discrete curl - not defined in 1D
% I = Nx by Nx identity matrix
% M = Nx+1 by Nx discrete mean matrix
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> [D,G,C,I,M]=build_ops(Grid);

Nx = Grid.Nx; Nfx = Grid.Nfx;

%% 1) Build sparse Divergence operator
B1 = ones(Nfx, 1);
B2 = -1.*B1; 
d = [0 1];
D = -1* 1/Grid.dx .* spdiags([B1 B2],d,Nx,Nfx);
full(D);

%% 2) Build sparse Gradient operator
%  Interior
B2 = -1.*B1; 
d = [0 -1];
G = 1/Grid.dx .* spdiags([B1 B2],d,Nfx,Nx);
% Set natural (homogeneous Neumann) boundary conditions
dof_f_bnd = [1 Nfx]; % all dof's on boundary
G(dof_f_bnd,:) = 0;
full(G);

% 3) Discrete Curl operator (not defined in 1D)
C = [];

% 4) Sparse Identity 
I = speye(Nx);
full(I);

% 5) Sparse Mean
% Interior
M = 0.5 .* spdiags([B1 B1],0,Nfx,Nx);
full(M);

% Boundaries
M(1,1) = 1; 
M(Nfx,Nx) = 1;
full(M);