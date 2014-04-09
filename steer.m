function su = steer(pts,degree)
%STEER create steering polynomials for HDTV penalties
%degree = degree of the HDTV penalty
%[PX,PY,PZ] = samples of sphere

PX = pts(1,:); PY = pts(2,:); PZ = pts(3,:);

if(degree == 1)
    su = [PX; PY; PZ];
elseif(degree == 2)
    su = [PX.*PX; PY.*PY; PZ.*PZ; 2*PX.*PY; 2*PX.*PZ; 2*PY.*PZ];
elseif(degree == 3)
    su = [PX.^3; PY.^3; PZ.^3; 3*(PX.^2).*PY; 3*(PX.^2).*PZ; 3*PX.*(PY.^2); 3*PZ.*(PY.^2); 3*PX.*(PZ.^2); 3*PY.*(PZ.^2); 6*PX.*PY.*PZ];
end

end