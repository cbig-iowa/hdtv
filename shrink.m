function Y = shrink(X,tau)
%Returns thresholded magnitudes of X.
Z = abs(X);
Z(Z==0) = 1;
Y = max(Z - tau, 0)./Z;

end


