clear all
close all
clc

%% Simple 2D Cartesian grids
% Cartesian grid
clf;
nx = 10;
ny = 10;
nz = 10;
G = cartGrid([nx, ny, nz], [nx ny nz]);
plotGrid(G);
view(3), colormap(jet), axis equal tight;
%% Heterogeneous permeability
rock                = makeRock(G, 0.5*milli*darcy, 0.2);
rock.perm(1:100)    = 1*milli*darcy();
rock.perm(300:400)  = 2*milli*darcy();
rock.perm(601:700)  = 3*milli*darcy();
rock.perm(901:1000) = 4*milli*darcy();
f(1)=figure(1);
plotCellData(G, log(rock.perm));
view(3), colormap(jet), axis equal tight; colorbar

%% Compute Geometry 
% Define faces, nodes, neighbours 
G          = computeGeometry(G);
%% Define Fluid

fluid      = initSingleFluid('mu' ,    1*centi*poise, ...
    'rho', 1014*kilogram/meter^3);
gravity reset on
%% Boundary conditions
bc  = pside([], G, 'TOP', 1000.*barsa());

%% Assemble and solve the linear system
% To solve the flow problem, we use the standard two-point
% flux-approximation method (TPFA), which for a Cartesian grid is the same
% as a classical seven-point finite-difference scheme for Poisson's
% equation. This is done in two steps: first we compute the
% transmissibilities and then we assemble and solve the corresponding
% discrete system.
T   = simpleComputeTrans(G, rock);
sol = simpleIncompTPFA(initResSol(G, 0.0), G, T, fluid, 'bc', bc);

%% Plot the face pressures
newplot;
plotFaces(G, 1:G.faces.num, convertTo(sol.facePressure, barsa()));
set(gca, 'ZDir', 'reverse'), title('Pressure [bar]')
view(3), colorbar


