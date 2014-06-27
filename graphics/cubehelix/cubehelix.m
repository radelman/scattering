function [map,lo,hi] = cubehelix(N,start,rots,hue,gamma,rng)
% Generate colormaps of Dave Green's Cubehelix colorscheme, with brightness control.
%
% (c) 2014 Stephen Cobeldick
%
% ### Function ###
%
% Returns a colormap with colors defined by Dave Green's Cubehelix colorscheme.
% The colormap nodes are selected along a tapered helix in the RGB color cube,
% with a continuous increase in perceived intensity. Black-and-white printing
% using postscript results in a monotonically increasing grayscale colorscheme.
%
% The scheme is defined here: http://astron-soc.in/bulletin/11June/289392011.pdf
% For more information and examples: http://www.mrao.cam.ac.uk/~dag/CUBEHELIX/
%
% The function has an optional input to control the colormap brightness range:
% this is useful for printing or matching an existing document colorscheme.
%
% Syntax:
%  map = cubehelix;
%  map = cubehelix(N);
%  map = cubehelix(N,start,rots,hue,gamma);
%  map = cubehelix(N,start,rots,hue,gamma,rng);
%  map = cubehelix([],...)
% [map,lo,hi] = cubehelix(...)
%
% See also CUBEHELIX_VIEW CUBEHELIX_FIND BREWERMAP RGBPLOT3 RGBPLOT COLORMAP COLORBAR SURF CONTOURF IMAGE CONTOURCMAP JET LBMAP
%
% ### Examples ###
%
% % Plot a scheme's RGB values:
% rgbplot(cubehelix(100))
%
% % Plot a scheme in an RGB cube:
% X = cubehelix(100,3,1.5,1,1);
% axes('ColorOrder',X, 'NextPlot','replacechildren', 'View',[-40,40])
% plot3(X(:,[1,1]).',X(:,[2,2]).',X(:,[3,3]).','.','MarkerSize',36)
% grid on
%
% % New colors for the "colormap" example:
% load spine
% image(X)
% colormap(cubehelix)
%
% % New colors for the "surf" example:
% [X,Y,Z] = peaks(30);
% surfc(X,Y,Z)
% colormap(cubehelix([],0.5,-1.5,1,1,[0.29,0.92]))
% axis([-3,3,-3,3,-10,5])
%
% ### Input and Output Arguments ###
%
% Inputs (*=default):
%  N     = NumericScalar, N>=0, an integer to define the colormap length.
%        = *[], uses the length of the current figure's colormap.
%  start = NumericScalar, *0.5,  the start color, with R=1, G=2, B=3 etc. (modulus 3).
%  rots  = NumericScalar, *-1.5, the number of R->G->B rotations over the scheme length.
%  hue   = NumericScalar, *1,    controls how saturated the colors are.
%  gamma = NumericScalar, *1,    can be used to emphasize low or high intensity values.
%  rng   = NumericVector, *[0,1], brightness levels of the colormap's endnodes. Size 1x2.
%
% Outputs:
%  map = NumericMatrix, size N-by-3, a colormap of RGB values between 0 and 1.
%  lo  = LogicalMatrix, same size as <map>, true where RGB values<0 were clipped.
%  hi  = LogicalMatrix, same size as <map>, true where RGB values>1 were clipped.
%
% [map,lo,hi] = cubehelix(N,start,rots,hue,gamma,rng)

% ### Input Checking ###
%
if nargin==0 || isempty(N)
    N = size(get(gcf,'colormap'),1);
else
    assert(isscalar(N)&&isreal(N),'Input <N> must be a real scalar numeric.')
end
if nargin<=1
    % Default parameters.
    start = 0.5;
    rots  = -1.5;
    hue   = 1;
    gamma = 1;
elseif nargin<=3
    % Parameters in a vector.
    if nargin==3
        rng = rots;
    end
    assert(numel(start)==4&&isreal(start),'Parameters must be a 1x4 real numeric vector.')
    gamma = start(4);
    hue =   start(3);
    rots =  start(2);
    start = start(1);
else
    % Parameters as individual scalar values.
    assert(isscalar(start)&&isreal(start),'Input <start> must be a real scalar numeric.')
    assert(isscalar(rots) &&isreal(rots), 'Input <rots> must be a real scalar numeric.')
    assert(isscalar(hue)  &&isreal(hue),  'Input <hue> must be a real scalar numeric.')
    assert(isscalar(gamma)&&isreal(gamma),'Input <gamma> must be a real scalar numeric.')
end
if nargin<=5&&nargin~=3
    rng = [0,1];
else
    assert(numel(rng)==2&&isreal(rng),'Input <rng> must be a 1x2 real numeric vector.')
end
%
% ### Core Function ###
%
map = zeros(N,3);
cof = [-0.14861,1.78277;-0.29227,-0.90649;1.97294,0];
%
for n = 1:N
    ang = 2*pi * (start/3+1+rots*(n-1)/(N-1));
    fra = ((n-1)/(N-1))^gamma;
    amp = hue * fra * (1-fra) / 2;
    fra = rng(1)*(1-fra) + rng(2)*fra;
    map(n,:) = fra + amp * (cof*[cos(ang);sin(ang)]);
end
%
lo = map<0;
hi = map>1;
map = max(0,min(1,map));
%
end
%----------------------------------------------------------------------END:cubehelix