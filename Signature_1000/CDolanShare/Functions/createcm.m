function [ cmap , Y ] = createcm( V, N )
%CREATECM Interpolates a colormap from colors given in V
%   cmap = CREATECM(V) creates a 64-by-3 matrix that corresponds to an
%   interpolated colormap from the colors in V, must be an RGB array of at
%   least 2 colors (where values range from 0 to 1).
%
%   [cmap,Y] = CREATECM(V) returns the relative luminance, Y, which is the
%   perceived lightness of a color. This is useful for checking how a
%   colormap might appear in greyscale.
%
%   cmap = CREATECM(V,N) creates a N-by-3 matrix, rather than the default
%   64-by-3 matrix for the colormap.
%
%   Examples:
%       V = [10 10 30; 20 50 75; 95 95 95; 75 30 20; 30 0 10]/100
%       cmap = createcm(V)          % returns a blue-white-red colormap
%
%       V = [0 15 20; 35 25 60; 75 40 50; 100 60 25; 95 100 35]/100
%       [cmap,Y] = createcm(V,128)  % returns a blue-purple-yellow colormap
%       contourf(randn(10))         % plot random colors
%       colormap(cmap)              % changes colormap to custom color
%       plot(Y)                     % shows relative luminance
%
%       V = [15 10 45; 20 60 55; 100 95 60]/100 % blue-green-yellow
%
%   A.P. Garcia

%%

if nargin == 1
    N = 64;     % default size for colormaps
end

m = size(V,1);
cmap = interp2(V, repmat(1:3,N,1), repmat(1:(m-1)/(N-1):m,3,1)');
Y = 0.2126*cmap(:,1)+0.7152*cmap(:,2)+0.0722*cmap(:,3);

end

