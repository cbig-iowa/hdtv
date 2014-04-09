function P = phantom3d_linear( N )
%Piecewise linear 3D Shepp-Logan

I = phantom3d(N);
[X,Y,Z] = meshgrid(1:N,1:N,1:N);
L = 2*X+Y+Z;
L = L/max(L(:));
P = L.*I; %multiply by simple linear function for a linear fade.

end

