% Nonlocal Function G_1 and G_2
function out = Gr(x,r)
    out = x.*(1-abs(x)/r).*(abs(x)<=r); 
    % loc = x.*(abs(x)<=r);
end