function [pts,weights] = getLebedevQuad(degree)
    %allowed = [ 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810 ];    

    leb = getLebedevSphere(degree);
    
    x = leb.x;
    y = leb.y;
    z = leb.z;
    w = leb.w;
    
    %take only samples in half-sphere:
    mask = (x > 0) | ((x==0) & (y>0)) | (z == 1);

    x = x(mask);
    y = y(mask);
    z = z(mask);
    w = w(mask);
        
    pts = [x'; y'; z'];
    weights = w';
end