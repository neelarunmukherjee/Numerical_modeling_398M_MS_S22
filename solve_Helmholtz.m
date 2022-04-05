function [hD,pD,qD] = solve_Helmholtz(D,G,I,M,phi,n,m,Grid,B,N,fn,BC,Zc,Gamma) % class
% author: you
% date: today

% Input: (to be completed)
% D = N by Nf discrete divergence operator
% G = Nf by N discrete gradient operator
% I = N by N identity
% phi = Ny by Nx matrix with porosity field
% n = permeability exponent
% m = compaction exponent
% B = Nc by N constraint matrix
% N = N by N-Nc basis for nullspace of constraint matrix
% fn = r.h.s. vector for Neuman BC
% BC = struture with info for boundary condtions
% Zc = Ny by Nx matrix containing the vertical coordinate of cell centers
% Gamma = N by 1 vector containing melting term

% Output:
% hD = N by 1 column vector of dimensionless overpressure heads
% pD = N by 1 column vector of dimensionless overpressure
% qD = Nf by 1 column vector of dimensionless relative fluid fluxes

% Porosity Matrices
Phi_n = phi.^n;
phi_n =comp_mean(Phi_n,M,1,Grid,1);
Phi_m = phi.^m;
phi_m =spdiags(Phi_m,0,Grid.Nx,Grid.Nx);

%% Solve mod. Helmholtz equations 

L  = -D*phi_n*G + phi_m;      % system matrix
fs = phi_m*Zc + Gamma;     % r.h.s.

%% Solve boundary value problem
hD = solve_lbvp(L,fs+fn,B,BC.g,N);
flux = @(h) -phi_n*G*h;
res = @(h,cell) L(cell,:)*h - fs(cell);
qD = comp_flux_gen(flux,res,hD,Grid,BC);
pD = hD-Zc;
end