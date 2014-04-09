function XYZ = tensor_prod( x,y,z )
%TENSOR_PROD 3d tensor product of vectors x,y,z.

sx = length(x);
sy = length(y);
sz = length(z);

x = reshape(x,[1,sx,1]);
y = reshape(y,[sy,1,1]);
z = reshape(z,[1,1,sz]);

X = repmat(x,[sy,1,sx]);
Y = repmat(y,[1,sx,sz]);
Z = repmat(z,[sy,sx,1]);

XYZ = X.*Y.*Z;

end

