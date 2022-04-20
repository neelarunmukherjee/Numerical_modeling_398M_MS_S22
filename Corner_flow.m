clear; close all; 
mu = 1;

%% Build staggered grids
Gridp.xmin = 0; Gridp.xmax = 1; Gridp.Nx = 100;
Gridp.ymin = 0; Gridp.ymax = 1; Gridp.Ny = 100;
Grid = build_stokes_grid(Gridp);

%% Build Stokes operators

[D,Edot,Dp,Gp,Z,I] = build_stokes_ops(Grid);
A = 2*mu*D*Edot; 
L = [A, -Gp;...
     Dp, Z];
fs = spalloc(Grid.N,1,0);

%% Build BC's

BC.dof_dir = [Grid.dof_ymax_vt;... % tangential velocity on the top
                Grid.dof_xmin_vt(2:end-1);... %no slip xmin
                Grid.dof_xmax_vt(2:end-1);... %no slip xmax
                Grid.dof_ymin_vt(2:end-1); %no slip ymin
                Grid.dof_pene;...  %no penetration
                Grid.dof_pc];      %pressure constraint
            
BC.g       = [ones(Grid.p.Nx-1,1);... % tangential velocity on the top
                zeros(length(Grid.dof_xmax_vt)-2,1);... %no slip xmin
                zeros(length(Grid.dof_xmax_vt)-2,1);... %no slip xmax
                zeros(length(Grid.dof_xmax_vt)-2,1);... %no slip ymin
                zeros(Grid.N_pene,1);... %no penetration
                0]; %pressure
[B,N,fn] = build_bnd(BC,Grid,I);

%% Solve for Stokes flow
u = solve_lbvp(L,fs+fn,B,BC.g,N);
v = u(1:Grid.p.Nf); p = u(Grid.p.Nf+1:end);
PSI = comp_streamfun(v,Grid.p);

%% Plot solution

[Xp,Yp] = meshgrid(Grid.x.xc,Grid.y.yc);
psi_max = max(PSI(:));
contour(Xp,Yp,PSI,20,'k'), hold on
contour(Xp,Yp,PSI,psi_max*[.2:.2:1],'r-')
legend('main circulation','corner eddies','location','south')
axis square
xlabel('x','fontsize',14)
ylabel('z','fontsize',14)
set(gca,'xtick',[0:.2:1],'ytick',[0:.2:1])
xlim([0 1]), ylim([0 1])
ta = annotation('textarrow');
s = ta.FontSize;
ta.FontSize = 12;
ta.Position = [0.4000 .9500 0.2000 0.0000];
ta.VerticalAlignment = 'bot';
ta.HorizontalAlignment = 'right';
text(.45,1.06,'v_{lid}','fontsize',12)
