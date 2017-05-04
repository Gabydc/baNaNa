close all
clear all
clc

%% Set up grid and petrophysical data
% We use a Cartesian grid of size nx-by-ny with homogeneous petrophysical
% data: permeability layers of 1 mD and 10^-1 mD and porosity of 0.2.

% Define the domain
pw = 1;
sz = 32*2^pw;
nxi = 1; nyi = 1; nx = sz; ny = sz;
Lx = sz; Ly = sz;
% Exponent of the permeability
per = 1;
% Values for the linear solver, tolerance, number of iterations
tol = 10^-5;
iter = 500;

%% Set up Grid
G = cartGrid([sz, sz, 1], [64, 64, 1]);
G = computeGeometry(G);
% Disable gravity
gravity off

%% Set up layers of permeability and constant porosity
rock.perm = ones(G.cells.num, 1)*1*milli*darcy;
%inhomogeneus permeability
lsize = round(sz*sz/8);
rock.poro = ones(G.cells.num, 1)*0.2;
for i = 1 : 2 : 8
    rock.perm(1+lsize*(i-1):lsize*i)  = repmat(10^(-per)*milli*darcy(), [lsize, 1]);
end


f(1) = figure(1);
clf
plotCellData(G, log(rock.perm));
view(0,90), colormap(jet), axis equal tight;


%% Compute half transmissibilities
% All we need to know to develop the spatial discretization is the reservoir
% geometry and the petrophysical properties. This means that we can compute
% the half transmissibilities without knowing any details about the fluid
% properties and the boundary conditions and/or sources/sinks that will
% drive the global flow:
hT = simpleComputeTrans(G, rock);

%% Fluid model
% When gravity forces are absent, the only fluid property we need in the
% incompressible, single-phase flow equation is the viscosity. However, the
% flow solver is written for general incompressible flow and requires the
% evaluation of a fluid object that can be expanded to represent more
% advanced fluid models. Here, however, we only use a simple fluid object
% that requires a viscosity and a density (the latter is needed when gravity
% is present)3 bar
gravity reset off
fluid = initSingleFluid('mu' , 1*centi*poise, ...
    'rho', 1014*kilogram/meter^3);


%% Define wells properties, 5 wells

for s = 5
    well(1:4) = -1;
    well(5)   = 3;
    xi = rand(nx*ny,1);

    num = num2str(s);
    
    wtype    = {'bhp', 'bhp', 'bhp', 'bhp', 'bhp'};
    wtarget  = [well(1),   well(2),   well(3),   well(4), well(5)] .* barsa();
    wrad     = [0.125, 0.125, 0.125, 0.125, 0.125] .* meter;
    wloc     = [  nxi,   nxi,     nx,   nx, nx/2;
        nyi,   ny,     nyi,   ny, ny/2];
    wname    = {'W1', 'W2', 'W3', 'W4', 'W5'};
    sgn      = [ 1 ,  1 ,  1 ,  1 ,1 ];
    W = [];
    for w = 1 : numel(wtype),
        W = verticalWell(W, G, rock, wloc(1,w), wloc(2,w), [], ...
            'Type', wtype{w}, 'Val', wtarget(w), ...
            'Radius', wrad(w), 'Name', wname{w}, ...
            'Sign', sgn(w), 'InnerProduct', 'ip_tpf');
    end
    figure(1)
    for i=1:numel(W)
 plotWell(G, W(i));
    end
view(0,90), camproj perspective, axis tight off
    
    %% Initialize state
    sol = initState(G, W, 0);
    
    
    %% Select solvers
    mrstModule add agmg
    lsolver = 3;
    switch lsolver
        case 1
            solver = AGMGSolverAD('tolerance', tol,'maxIterations', iter);
            ls = 'AGMG';
        case 2
            solver = GMRES_ILUSolverAD('tolerance', tol,'maxIterations', iter);
            ls = 'GMRES';
        case 3
            solver = PCG_ICSolverAD('tolerance', tol,'maxIterations', iter);
            ls = 'PICCG';
        case 4
            solver = DPCG_ICSolverAD('tolerance', tol,'maxIterations', iter);
            ls = 'DPICCG';
        case 5
            solver = BackslashSolverAD();
            ls = 'Backslash';
    end
    
    fn = @(A, b) solver.solveLinearSystem(A, b);
    psolve = @(state) incompTPFA_g_o(state, G, hT, fluid, 'wells', W,'MatrixOutput',true,'LinSolve', fn,'verbose',true);
    [sol,report] = psolve(sol);
    p = sol.pressure;
    
    A = sol.A(1:G.cells.num,1:G.cells.num);
    b = sol.rhs(1:G.cells.num);
    xb = A\b;
    
    
    figure(1)
    [ht] = plotingsolution(G,W,[ls], p,1) ;
    colorbar
    [ht] = plotingsolution(G,W,'Backslash',xb,2);
    colorbar
      return
    
    %% Use matrix

    
    l=ichol(A);
    [x11,hc0,hc01,hc02]=ICCG_0(A,b,xi,iter,tol,l,0);
    figure
    plotGrid(G, 'FaceAlpha', 0, 'EdgeAlpha', .1);
 plotCellData(G,x11)
view(0,90)

axis equal tight; colormap(jet(128));
title('ICCG');

  
    
   
end



