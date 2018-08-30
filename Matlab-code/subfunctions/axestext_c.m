function handle=axestext_c(xpos,ypos,string,varargin)

% axestext(xpos,ypos,string,varargin)
% James R code
% Write a text string in a specified location relative to current axes.
%  Required Inputs:
%   xpos, ypos  - where, relative to axes ((0,0)=BL, (1,1)=TR, etc)
%   string      - The input string
%  Optional (VARARGIN) Inputs:
%	any valid arguments for 'text', PLUS 'blank' for blanking under tex

if (nargin<3), error('Must supply xpos,ypos,string'), end

xl=xlim; xd=diff(xl); xinc=xd./100;
yl=ylim; yd=diff(yl); yinc=yd./100;

haseq={'right';'left';'center';'right';'left'};
vaseq={'top';'bottom';'middle';'top';'bottom'};
amult=[-1 1 0 -1 1];
xarea=bin_c(xpos,[0 1/4 3/4 1]); yarea=bin_c(ypos,[0 1/4 3/4 1]);

xp=xl(1)+xpos*xd+amult(xarea)*xinc;
yp=yl(1)+ypos*yd+amult(yarea)*yinc;

halign=haseq{xarea};
valign=vaseq{yarea};

k=find(strncmp(varargin,'horiz',5));
if isempty(k), varargin={varargin{:},'horizontal',halign}; end
k=find(strncmp(varargin,'vert',4));
if isempty(k), varargin={varargin{:},'vertical',valign}; end
h=puttext(xp,yp,string,varargin{:});
if nargout>0, handle=h;, end

