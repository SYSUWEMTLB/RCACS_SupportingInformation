function quiver2(x,y,s,xcoord,ycoord)

% QUIVER2 Improved quiver to allow x,y specifications: superceded by quiver_ext
% 
%QUIVER2   1.1   92/04/14
%	Modification to the standard QUIVER plot to allow specification
%	of the x and y coordinates -- this enables a quiver plot to
%	be drawn with meaningful axes, and in particluar to be overlaid
%	(using HOLD ON) upon a CONTOUR plot.
%
%	QUIVER2(Z) draws a graph that displays the angle and magnitude
%	of the complex elements of Z as arrows emanating from equally
%	spaced grid positions, one arrow per matrix element.
%
%	QUIVER2(X,Y) is equivalent to QUIVER2(X+i*Y).  It displays the
%	quiver plot for the angles and magnitudes of the elements of
%	matrices X and Y.
%
%	QUIVER2(Z,'S') and QUIVER2(X,Y,'S') use line style 'S' where
%	'S' is any legal linetype as described under the PLOT command.
%
%	QUIVER2(Z,'S',XCOORD,YCOORD) and QUIVER2(X,Y,'S',XCOORD,YCOORD) 
%	use XCOORD and YCOORD as the axis scales.  For example, vectors
%	giving the longitude and latitude coordinates respectively of
%	points at which the X,Y data are defined
%
%	See also QUIVER, GRADIENT, COMPASS, FEATHER, and ROSE.
%
%	BUGS: Doesn't seem to like plotting subsets, 
%	e.g. QUIVER2(U(Y,X),V(Y,X),XCOORD(X),YCOORD(Y))?
%
%	John Wilkin 24/3/92
%	adapted from the MATLAB QUIVER function

%	Charles R. Denham, MathWorks 3-20-89
%	Copyright (c) 1989 by the MathWorks, Inc.

xx = [0 1 .8 1 .8].';
yy = [0 0 .08 0 -.08].';
arrow = xx + yy.*sqrt(-1);

if nargin == 1,
   s = 'r-';
   y = imag(x); x = real(x);
   [m,n] = size(x);
   xcoord = 1:n;
   ycoord = 1:m;
elseif nargin == 2,
   if isstr(y)
      s = y;
      y = imag(x); x = real(x);
   else
      s = 'r-';
   end
   [m,n] = size(x);
   xcoord = 1:n;
   ycoord = 1:m;
elseif nargin == 3
   if isstr(s)
      [m,n] = size(x);
      xcoord = 1:n;
      ycoord = 1:m;
   else
      ycoord = s;
      xcoord = y;
      y = imag(x); x = real(x);
      s = 'r-';
   end
elseif nargin == 4
   if isstr(y)
      ycoord = xcoord;
      xcoord = s;
      s = y;
      y = imag(x); x = real(x);
   else
      ycoord = xcoord;
      xcoord = s;
      s = 'r-';
   end
end

[xx,yy] = meshdom(xcoord, ycoord);
grid = xx + yy.*sqrt(-1); grid = grid(:);
x = x(:); y = y(:);
z = (x + y.*sqrt(-1)).';
scale = 0.90 ./ max(max(abs(z)));
a = scale * arrow * z + ones(5,1) * grid.';
plot(real(a), imag(a), s);
