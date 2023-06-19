%{

The methods of inpaint_nans

Digital inpainting is the craft of replacing missing elements in an
"image" array. A Google search on the words "digita inpainting" will turn 
up many hits. I just tried this search and found 18300 hits.

If you wish to do inpainting in matlab, one place to start is with my
inpaint_nans code. Inpaint_nans is on the file exchange:

http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=4551&objectType=file

It looks for NaN elements in an array (or vector) and attempts to interpolate
(or extrapolate) smoothly to replace those elements.

The name "inpainting" itself comes from the world of art restoration.
Damaged paintings are restored by an artist/craftsman skilled in matching
the style of the original artist to fill in any holes in the painting.

In digital inpainting, the goal is to interpolate in from the boundaries
of a hole to smoothly replace an artifact. Obviously, where the hole is
large the digitally inpainted repair may not be an accurate approximation
to the original.

Inpaint_nans itself is really only a boundary value solver. The basic idea
is to formulate a partial differential equation (PDE) that is assumed to
apply in the domain of the artifact to be inpainted. The perimeter of the
hole supplies boundary values for the PDE. Then the PDE is approximated
using finite difference methods (the array elements are assumed to be
equally spaced in each dimension) and then a large (and very sparse) linear
system of equations is solved for the NaN elements in the array.

I've chosen a variety of simple differental equation models the user can
specify to be solved. All the methods current use a basically elliptic
PDE. This means that the resulting linear system will generally be well
conditioned. It does mean that the solution will generally be fairly smooth,
and over large holes, it will tend towards an average of the boundary
elements. These are characteristics of the elliptic PDEs chosen. (My hope
is to expand these options in the future.)

%}

%%

% Lets formulate a simple problem, and see how we could solve it using
% some of these ideas.
A = [0 0 0 0;1 NaN NaN 4;2 3 5 8];

% Although we can't plot this matrix using the functions surf or mesh,
% surely we can visualize what the fudamental shape is.

% There are only two unknown elements, the artifacts that inpaint_nans
% would fill in: A(2,2) and A(2,3).

% For an equally spaced grid, the Laplacian equation (or Poisson's equation
% of heat conduction at steady state if you prefer. Or, for the fickle,
% Ficke's law of diffusion would apply.) All of these result in the PDE
%
%   u_xx + u_yy = 0
%
% where u_xx is the second partial derivative of u with respect to x,
% and u_yy is the second partial with respect to y. 
%
% Approximating this PDE using finite differences for the partial
% derivatives, implies that at any node in the grid, we could replace
% it by the average of its 4 neighbors. Thus the two NaN elements
% generate two linear equations:
%
%  A(2,2) = (A(1,2) + A(3,2) + A(2,1) + A(2,3)) / 4
%  A(2,3) = (A(1,3) + A(3,3) + A(2,2) + A(2,4)) / 4
%
% Since we know all the parameters but A(2,2) and A(2,3), substitute their
% known values.
%
%  A(2,2) = (0 + 3 + 1 + A(2,3)) / 4
%  A(2,3) = (0 + 5 + A(2,2) + 4) / 4
%
% Or,
%
%  4*A(2,2)  - A(2,3) = 4
%  -A(2,2) + 4*A(2,3) = 9
%
% We can solve for the unkowns now using 
u = [4 -1;-1 4]\[4;9]

A(2,2) = u(1);
A(2,3) = u(2);

% and finally plot the surface
close
surf(A)
title 'A simply inpainted surface'

% Neat huh? For an arbitrary number of NaN elements in an array,
% the above scheme is all there is to method 2 of inpaint_nans,
% together with a very slick application of sparse linear algebra
% in Matlab.

% Method 0 is very similar, but I've optimized it to build as
% small a linear system as possible for those cases where an array
% has only a few NaN elements.

% Method 1 is another subtle variation on this scheme, but it
% tries to be slightly smoother at some cost of efficiency, while
% still not modifying the known (non-NaN) elements of the array.

% Method 5 of inpaint_nans is also very similar to method 2, except
% that it uses a simple average of all 8 neighbors of an element.
% Its not actually an approximation to our PDE.

% Method 3 is yet another variation on this theme, except the PDE
% model used is one more suited to a model of a thin plate than for
% heat diffusion. Here the governing PDE is:
%
%   u_xxxx + 2*u_xxyy + u_yyyy = 0
%
% again discretized into a linear system of equations. 

%%

% Finally, method 4 of inpaint_nans has a different underlying
% model. Pretend that each element in the array was connected to
% its immediate neighbors to the left, right, up, and down by
% "springs". They are also connected to their neighbors at 45
% degree angles by springs with a weaker spring constant. Since
% the potential energy stored in a spring is proportional to its
% extension, we can formulate this again as a linear system of
% equations to be solved. For the example above, we would generate
% the set of equations:

%  A(2,2) - A(1,2) = 0
%  A(2,2) - A(2,1) = 0
%  A(2,2) - A(3,2) = 0
%  A(2,2) - A(2,3) = 0
% (A(2,2) - A(1,1))/sqrt(2) = 0
% (A(2,2) - A(1,3))/sqrt(2) = 0
% (A(2,2) - A(3,1))/sqrt(2) = 0
% (A(2,2) - A(3,3))/sqrt(2) = 0
%  A(2,3) - A(1,3) = 0
%  A(2,3) - A(2,2) = 0
%  A(2,3) - A(3,3) = 0
%  A(2,3) - A(2,4) = 0
% (A(2,3) - A(1,2))/sqrt(2) = 0
% (A(2,3) - A(1,4))/sqrt(2) = 0
% (A(2,3) - A(3,2))/sqrt(2) = 0
% (A(2,3) - A(3,4))/sqrt(2) = 0

% Substitute for the known elements to get

%  A(2,2) - 0 = 0
%  A(2,2) - 1 = 0
%  A(2,2) - 3 = 0
%  A(2,2) - A(2,3) = 0
% (A(2,2) - 0)/sqrt(2) = 0
% (A(2,2) - 0)/sqrt(2) = 0
% (A(2,2) - 2)/sqrt(2) = 0
% (A(2,2) - 5)/sqrt(2) = 0
%  A(2,3) - 0 = 0
%  A(2,3) - A(2,2) = 0
%  A(2,3) - 5 = 0
%  A(2,3) - 4 = 0
% (A(2,3) - 0)/sqrt(2) = 0
% (A(2,3) - 0)/sqrt(2) = 0
% (A(2,3) - 3)/sqrt(2) = 0
% (A(2,3) - 8)/sqrt(2) = 0

% This system is also solvable now:
r2 = 1/sqrt(2);
M=[1 0;1 0;1 0;1 -1;r2 0;r2 0;r2 0;r2 0;0 1;-1 1;0 1;0 1;0 r2;0 r2;0 r2;0 r2];
v = M\[0 1 3 0 0 0 2*r2 5*r2 0 0 5 4 0 0 3*r2 8*r2]'

A(2,2) = v(1);
A(2,3) = v(2);

% and finally plot the surface
surf(A)
title 'A simply inpainted surface using a spring model'

%%

% Why did I provide this approach, based on a spring metaphor?
% As you should have observed, methods 2 and 4 are really quite close
% in what they do for internal NaN elements. Its on the perimeter that
% they differ significantly. The diffusion/Laplacian model will
% extrapolate smoothly, and as linearly as possible. The spring model
% will tend to extrapolate as a constant function.