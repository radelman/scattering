function [vec,rng,resnorm] = cubehelix_find(map)
% Retrieve the parameter values from an existing Cubehelix colorscheme colormap.
%
% (c) 2014 Stephen Cobeldick
%
% So you have a nice Cubehelix colormap, but you can't remember the exact
% parameter values that were used to define it... this function can help!
% Given an RGB colormap, uses a least-squares solver to approximate the original
% Cubehelix colorschemes' parameter values, and returns these in a vector.
%
% Syntax:
%  vec = cubehelix_find(map)
%  [vec,rng] = cubehelix_find(map)
%
% Note1: Requires the Optimization Toolbox function "lsqnonlin".
% Note2: The parameter <start> is modulus three, i.e. 3==0.
%
% See also CUBEHELIX CUBEHELIX_VIEW BREWERMAP RGBPLOT3 RGBPLOT COLORMAP LSQNONLIN OPTIMSET
%
% ### Examples ###
%
% cubehelix_find(cubehelix(10))
%  ans = [0.5,-1.5,1,1]
%
% map = cubehelix(10,1.4,-0.7,0.9,1.2);
% cubehelix_find(map)
%  ans =  [1.4,-0.7,0.9,1.2]
%
% map = cubehelix(64,0.2,1.6,1,0.8,[0.29,0.92]);
% [vec,rng] = cubehelix_find(map)
%  vec = [0.2,1.6,1,0.8]
%  rng = [0.29,0.92]
% all(cubehelix(64,vec,rng)-map <= 1e-10)
%  ans = [true,true,true]
%
% ### Input and Output Arguments ###
%
% Inputs:
%  map = NumericMatrix, a colormap of a Cubehelix colorscheme.
% Outputs:
%  vec = NumericVector, [start,rots,hue,gamma], the parameter values derived from <map>.
%  rng = NumericVector, [lb,ub], brightness levels at <map>'s endnodes.
%
% [vec,rng] = cubehelix_find(map)

siz = size(map);
assert(isfloat(map)&&numel(siz)==2&&siz(2)==3&&all(0<=map(:)&map(:)<=1),...
    'Input argument <map> must be a colormap of RGB values (size Nx3).')
%
% Parameter bounds:
lb = [0,-20,0,0];
ub = [3,+20,3,3];
opt = optimset('Display','off');
%
% Grayscale colormap equivalent:
mag = sum(map*[0.298936;0.587043;0.114021],2);
rng = mag([1,end]).';
%
N = size(map,1);
% Estimate gamma from grayscale:
g = lsqnonlin(@(g)mag-(diff(rng)*((0:N-1).'/(N-1)).^g+rng(1)),1,[],[],opt);
% Estimate rotations from the number of peaks&troughs:
D = diff(bsxfun(@minus,map,mag),1,1);
idx = D(1:end-1,:)>=0 & D(2:end,:)<0;
idy = D(1:end-1,:)<0 & D(2:end,:)>=0;
rots = ceil(mean(sum([idx,idy],1)));
%
% Define the solver function:
F = @(p)map-cubehelix(N,p,rng);
% Range of start and rotation values to try:
[S,R] = meshgrid(0:(1/(1+rots)):2.99,[-rots,rots]);
% Solve!
[Z,rn] = arrayfun(@(s,r)lsqnonlin(F,[s,r,1,g],lb,ub,opt),S,R,'UniformOutput',false);
% Pick the best solution:
[~,idz] = min([rn{:}]);
resnorm = rn{idz};
vec = Z{idz};
% Start value is modulus three:
vec(1) = mod(vec(1),3);
%
end
%----------------------------------------------------------------------END:cubehelix_find