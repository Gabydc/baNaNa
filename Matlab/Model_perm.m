%% Simple 2D Cartesian grids
% Cartesian grid
clf;
nx = 10;
ny = 10;
nz = 10;
G = cartGrid([nx, ny, nz], [nx ny nz]);
plotGrid(G);
view(3), colormap(jet), axis equal tight;

%% Define the model
% To set up a model, we need: a grid, rock properties (permeability), a
% fluid object with density and viscosity, and boundary conditions.

G          = computeGeometry(G);
% Homogeneous permeability
rock = makeRock(G, 1*milli*darcy, 0.2);
f(1)=figure(1);
plotCellData(G, log(rock.perm));
view(3), colormap(jet), axis equal tight; colorbar
%% Heterogeneous permeability
rock.perm(301:400)  = 2*milli*darcy();
rock.perm(601:700)  = 3*milli*darcy();
f(1)=figure(1);
plotCellData(G, log(rock.perm));
view(3), colormap(jet), axis equal tight;
%%  Gaussian distribution
p = gaussianField(G.cartDims, [1 100000]);
rock = makeRock(G,p(:),0.2);

f(1)=figure(1);
plotCellData(G, (rock.perm));


%% Pre-defined Models (SPE10)
nx = 60;
ny = 220;
layers = 1:10;
cartDims = [  nx,  ny,   numel(layers)];

G = cartGrid([nx ny numel(layers)]);
G = computeGeometry(G);
rock = getSPE10rock(layers);
perm = rock.perm(:,1);

f(1)=figure(1);
plotCellData(G, log(perm),'linestyle','none');
view(3), colormap(jet), axis equal tight;

%% Define Fluid

fluid      = initSingleFluid('mu' ,    1*centi*poise, ...
    'rho', 1014*kilogram/meter^3);
bc  = pside([], G, 'TOP', 100.*barsa());

gravity reset on

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
set(gca,'DataAspect',[1 1 10]);
