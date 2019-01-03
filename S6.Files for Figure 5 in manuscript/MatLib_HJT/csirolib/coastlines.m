function coastlines(r,symbol)
% COASTLINES Draws coast lines
%
%       Revision: 1.8  Date: 94/05/03
%	function coastlines(r,symbol)
%
%	Draws coastlines.  It assumes that the x and y axes have been
%	set by an earlier call to contour or xcontour and refer to
%	longitude and latitude.
% 
%	COASTLINES(R) plots every R'th subsample of the coastline data. 
%	Default is R = 2.  Maximum resolution (R=1) has 48352 data pairs.
%
%	COASTLINES('SYMBOL') and COASTLINES(R,'SYMBOL') uses the string 
%	'SYMBOL' to plot the coastline points.  Default is 'b.'.  COASTLINES 
%	performs no "pen-up" moves while draws so a solid linetype 
%	(e.g. 'SYMBOL' = '-') may draw undesired lines across the plot. 
%
%	EXAMPLE:  Plot coastlines over a contour plot created with xcontour.
%	xcontour(z,v,long,lat,badflag); hold on; coastlines(1,'ro');
%
%	FILES: Requires the file coastlines.mat which contains the
%	Reid coastline data.
%
%       CALLER:   general purpose
%       CALLEE:   none
%
% 	AUTHOR: John Wilkin 23/3/92  & hacked a bit by Jim Mansbridge
%               4/6/92 and 13/12/93.
%=======================================================================

%       @(#)coastlines.m   1  1.3
%
%-----------------------------------------------------------------------

if nargin < 1			% defaults
  r = 2;
  symbol = 'b.';
end

if nargin == 1,
  if isstr(r) == 1,
    symbol = r;
    r = 2;                      % input is plot SYMBOL, set default R
  else
    symbol = 'b.';		% input is R, set default plot SYMBOL
  end
end  

subset = 1:r:48352;

load coastlines;
val_axis = axis;

if val_axis(1) < 0             % Plot coastlines with -180 < lon < 180
  xx = lon(subset,:);
  index = find(xx > 180);
  xx(index) = xx(index) - 360;
  plot(xx,lat(subset,:),symbol);
end

if val_axis(2) > 0           % Plot coastlines with 0 < lon < 360
  plot(lon(subset,:),lat(subset,:),symbol);
end
