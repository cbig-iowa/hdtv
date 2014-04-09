function out = BsplineDiff(x,order,degree)

if(degree > order)
    disp('degree of derivative cannot be greater than order');
    return;
end
shift = degree/2:-1:-degree/2;

filterinit = [1,-1];filter=1;
for i=1:degree,
    filter = conv(filter,filterinit);
end

out = zeros(size(x));
for i=1:length(filter),
    out = out + filter(i)*Bspline(x+shift(i),order-degree);
end
