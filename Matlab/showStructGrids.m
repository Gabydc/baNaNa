%% Simple 2D Cartesian grids
% Cartesian grid
clf;
nx = 10;
ny = 10;
nz = 10;
G = cartGrid([nx, ny, nz], [nx ny nz]);
plotGrid(G); 
view(3)
%% show a graded grid

clf;
dx=cos((-1:0.1:1)*pi); x =dx; 
%dx = 1-0.5*cos((-1:0.1:1)*pi); x = -1.15+0.1*cumsum(dx);
y = 0:0.05:1;
G = tensorGrid(x, sqrt(y));
plotGrid(G);


%% removeCells: Create an ellipsoid grid
clf;
x = linspace(-2,2,21);
G = tensorGrid(x,x,x);
subplot(1,2,1); plotGrid(G);view(3); axis equal

subplot(1,2,2); plotGrid(G,'FaceColor','none','LineWidth',0.01);
G = computeGeometry(G);
c = G.cells.centroids;
r = c(:,1).^2+ 0.5*c(:,2).^2;
%r = c(:,1).^2 + 0.25*c(:,2).^2+0.25*c(:,3).^2;
G = removeCells(G, find(r>1));
plotGrid(G); view(-70,70); axis equal

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
