function cubehelix_view(start,rots,hue,gamma,rng)
% Create an interactive figure for Cubehelix colormap parameter selection. With demo!
%
% (c) 2014 Stephen Cobeldick
%
% View any of Dave Green's Cubehelix colorschemes in a figure. A 2D lineplot
% and two colorbars show both the RGB colors and the equivalent grayscale colors.
%
% Syntax:
%  cubehelix_view
%  cubehelix_view(start,rots,hue,gamma)
%  cubehelix_view(start,rots,hue,gamma,rng)
%
% Six sliders allow real-time interactive adjustment of the Cubehelix
% parameter values, which are also displayed. The parameters can also be
% set/reset by calling the function with these values as input arguments.
%
% Warnings are displayed in the figure if any RGB values are clipped, or if
% the grayscale colormap is not strictly monontonic increasing/decreasing.
%
% Clicking the button labeled 'Demo' provides an endless display of randomly
% generated CubeHelix color schemes. Click again to stop this display.
%
% The scheme is defined here: http://astron-soc.in/bulletin/11June/289392011.pdf
% For more information and examples: http://www.mrao.cam.ac.uk/~dag/CUBEHELIX/
%
% See also CUBEHELIX CUBEHELIX_FIND BREWERMAP RGBPLOT3 RGBPLOT COLORMAP COLORMAPEDITOR COLORBAR UICONTROL ADDLISTENER
%
% ### Input Arguments ###
%
% Inputs (*=default):
%  start = NumericScalar, *0.5,  the start color, with R=1, G=2, B=3 etc. (modulus 3).
%  rots  = NumericScalar, *-1.5, the number of R->G->B rotations over the scheme length.
%  hue   = NumericScalar, *1,    controls how saturated the colors are.
%  gamma = NumericScalar, *1,    can be used to emphasize low or high intensity values.
%  rng   = NumericVector, *[0,1], brightness levels of the colormap's endnodes. Size 1x2.
%
% cubehelix_view(start,rots,hue,gamma,rng)

switch nargin
    case 0
        chvUpDt(0.5,-1.5,1,1,0,1);
    case 1
        assert(numel(start)==4,'First input can be a vector of four parameters.')
        chvUpDt(start(1),start(2),start(3),start(4),0,1);
    case 2
        assert(numel(start)==4,'First input can be a vector of four parameters.')
        assert(numel(rots)==2,'Second input can be a vector of two parameters.')
        chvUpDt(start(1),start(2),start(3),start(4),rots(1),rots(2));
    case 4
        chvUpDt(start,rots,hue,gamma,0,1);
    case 5
        assert(numel(rng)==2,'Fifth input can be a vector of two parameters.')
        chvUpDt(start,rots,hue,gamma,rng(1),rng(2));
    otherwise
        error('Wrong number of inputs. Enter parameters individually or in vectors.')
end
%
end
%----------------------------------------------------------------------END:cubehelix_view
function vals = chvUpDt(varargin)
% Draw a new figure or update an existing figure. Callback for sliders & demo.
%
persistent pltA pltL imgA imgI uicS txtS txtW
%
% LHS and RHS slider bounds/limits and step sizes:
lb = [0,-3, 0, 0, 0, 0];
rb = [3, 3, 3, 3, 1, 1];
sp = [1, 1, 1, 1, 1, 1;...
      5, 5, 5, 5, 2, 2].'./10;
%
switch nargin
    case 1 % Demo update
        vals = varargin{1};
    case 2 % Slider callback
        vals = arrayfun(@(h)get(h,'Value'),uicS);
    case 6 % Function call
        assert(all(cellfun(@isfloat,varargin)),'Parameters must be numeric values.')
        assert(all(cellfun(@isscalar,varargin)),'Parameters must be scalar values.')
        vals = [varargin{:}];
        % Create a new figure:
        if isempty(pltL) || ~all(ishghandle(pltL))
            [pltA,pltL,imgA,imgI,uicS,txtS,txtW] = chvPlot(lb,rb,sp);
        end
    otherwise
        error('This really should not happen...')
end
% Update slider positions:
if nargin~=2
    arrayfun(@(h,n,x,v)set(h,'Value',max(n,min(x,v))),uicS,lb,rb,vals);
end
%
% Update parameter value text:
arrayfun(@(h,v)set(h,'String',num2str(v,'%.2f')),txtS,vals);
% Update XData and axes limits:
N = 128+round(pow2(log2(1+abs(vals(2))),6));
set(pltA, 'XLim',[1,N]);
arrayfun(@(h)set(h, 'XData', 1:N), pltL);
arrayfun(@(h)set(h, 'YLim', [0,N]+0.5), imgA);
% Get Cubehelix colormap:
[map,lo,hi]  = cubehelix(N,vals(1:4),vals(5:6));
% Update images/colorbars and line data values:
mag = sum(map*[0.298936;0.587043;0.114021],2);
set(imgI(1), 'CData',reshape(map,N,[],3))
set(imgI(2), 'CData',repmat(mag,[1,1,3]))
set(pltL(1), 'YData',map(:,1))
set(pltL(2), 'YData',map(:,2))
set(pltL(3), 'YData',map(:,3))
set(pltL(4), 'YData',mag)
% Update warning text:
mad = diff(mag);
str = {'Not Monotonic';'Clipped'};
set(txtW,'String',str([any(mad<=0)&&any(0<=mad),any(lo(:))||any(hi(:))]));
%
end
%----------------------------------------------------------------------END:chvUpDt
function [pltA,pltL,imgA,imgI,uicS,txtS,txtW] = chvPlot(lb,rb,sp)
% Draw a new figure with RGBplot axes, ColorBar axes, and uicontrol sliders.
%
% Parameter names:
names = {'start';'rotations';'hue';'gamma';'range_L';'range_R'};
gap = 0.01;
hgt = 0.70;
lft = 0.20;
rgt = 0.24;
wdt = 1-lft-rgt-2*gap;
brh = (1-gap-hgt)/numel(names)-gap;
%
% RGBplot gray color:
gry = 0.7*[1,1,1];
% RGBplot X-values:
X = (1:2).'*[1,1,1,1];
%
% Add 2D lineplot:
figH = figure('HandleVisibility','callback', 'NumberTitle','off',...
    'Name','Cubehelix Interactive Parameter Selector', 'Color','white');
pltA = axes('Parent',figH, 'Position',[gap,1-hgt+gap,lft+wdt,hgt-2*gap], 'Visible','off',...
    'ColorOrder',[1,0,0;0,1,0;0,0,1;gry], 'XTick',[], 'YTick',[], 'YLim',[0,1]);
pltL = line(X,X, 'Parent',pltA);
% Add warning text:
txtW = text('Parent',pltA, 'Units','normalized', 'Position',[0,1],...
    'HorizontalAlignment','left', 'VerticalAlignment','top');
% Add demo button:
uicontrol(figH, 'Style','togglebutton', 'Units','normalized',...
    'Position',[wdt+lft/2+gap,1-hgt+gap,lft/2,brh], 'String','Demo',...
    'Max',1, 'Min',0, 'Callback',@chvDemo);
%
temp = [0,0,0,0,0,0];
txtS = [0,0,0,0,0,0];
uicS = [0,0,0,0,0,0];
for n = 1:6
    % Add parameter sliders:
    Y = gap+(6-n)*(brh+gap);
    uicS(n) = uicontrol(figH,'Style','slider', 'Units','normalized',...
        'Position',[lft+gap,Y,wdt,brh], 'Min',lb(n), 'Max',rb(n),...
        'SliderStep',sp(n,:)/(rb(n)-lb(n)));
    addlistener(uicS(n), 'Value', 'PostSet', @chvUpDt);
    % Add text to show slider parameter values:
    temp(n) = uicontrol(figH,'Style','text', 'Units','normalized',...
        'Position',[gap,Y,lft/2,brh],'String',names{n});
    txtS(n) = uicontrol(figH,'Style','text', 'Units','normalized',...
        'Position',[gap+lft/2,Y,lft/2,brh],'String','X');
end
%
% Add colorbars:
C = reshape([1,1,1],1,[],3);
imgA(1) = axes('Parent',figH, 'Visible','off', 'Units','normalized',...
    'Position',[1-rgt/1,gap,rgt/2-gap,1-2*gap], 'YLim',[0.5,1.5]);
imgA(2) = axes('Parent',figH, 'Visible','off', 'Units','normalized',...
    'Position',[1-rgt/2,gap,rgt/2-gap,1-2*gap], 'YLim',[0.5,1.5]);
imgI(1) = image('Parent',imgA(1), 'CData',C);
imgI(2) = image('Parent',imgA(2), 'CData',C);
%
end
%----------------------------------------------------------------------END:chvPlot
function chvDemo(tgh,~)
% While toggle button is depressed run a loop showing random Cubehelix schemes.
%
% Step size between updates:
step = 0.09;
% Initial values:
prms = chvUpDt([],[]);
ziel = prms;
% Functions to randomly find new values:
randfn(5:6) = {@()rand(1,1).^59,@()1-rand(1,1).^59};
randfn(3:4) = {@()sqrt(-log(rand(1,1))*2)};
randfn(1:2) = {@()3*rand(1,1),@()2*randn(1,1)};
%
while ishghandle(tgh)&&get(tgh,'Value')
    % While the toggle button is down, step values.
    isar = prms==ziel;
    ziel(isar) = round(100*cellfun(@(f)f(),randfn(isar)))/100;
    isnx = abs(ziel-prms)<=step;
    prms(isnx) = ziel(isnx);
    isgo = ~(isar|isnx);
    prms(isgo) = prms(isgo)+step*sign(ziel(isgo)-prms(isgo));
    % Update figure:
    prms = chvUpDt(prms);
    % Faster / slower ...?
    pause(0.11);
end
%
end
%----------------------------------------------------------------------END:chvDemo