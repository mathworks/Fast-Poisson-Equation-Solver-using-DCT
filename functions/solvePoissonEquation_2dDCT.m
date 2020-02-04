function [x,y,u] = solvePoissonEquation_2dDCT(Nx,Ny)
% Copyright 2020 The MathWorks, Inc.

% What is arguments?
% see: https://jp.mathworks.com/help/matlab/matlab_prog/argument-validation-functions.html
arguments 
    Nx (1,1) {mustBeInteger,mustBeGreaterThan(Nx,1)}
    Ny (1,1) {mustBeInteger,mustBeGreaterThan(Ny,1)}
end

dx = 2/Nx; dy = 2/Ny;
% The location of uÅiNOTE: Staggered GridÅj
xx = -1+dx/2:dx:1-dx/2;
yy = -1+dy/2:dy:1-dy/2;
[x,y] = ndgrid(xx,yy); % Note: not meshgrid

% modified wavenumber
kx = 0:Nx-1;
ky = 0:Ny-1;
mwx = 2*(cos(pi*kx/Nx)-1)/dx^2;
mwy = 2*(cos(pi*ky/Ny)-1)/dy^2;

% 2D DCT of f (Right hand side)
f = - 2*pi^2.*cos(pi*x).*cos(pi*y);
fhat = dct2(f); % Needs Image Processing Toolbox
% fhat = dct(dct(f)')'; % Same as above (Needs Signal Processing Toolbox instead)

[MWX, MWY] = ndgrid(mwx,mwy);
uhat = fhat./(MWX+MWY);

% The solution is not unique (uhat(0,0) = inf);
% Here we fix the mean ( with kx=0,ky=0) to be 0
uhat(1,1) = 0;

% Inverse 2D DCT
u = idct2(uhat); % Needs Image Processing Toolbox
% u = idct(idct(uhat)')'; % Same as above (Needs Signal Processing Toolbox instead)

end
