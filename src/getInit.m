% Define initial conditions
function init = getInit(x,p)
 
    noise = (rand(size(x))-0.5)*2; 
    noise = noise - mean(noise);
    uinit = ( ones(size(x))+0.3*noise ) * p.u0;
    % uinit = u0 + u0*0.1*sin(x*2*pi/L); 
    
    noise = (rand(size(x))-0.5)*2; 
    noise = noise - mean(noise);
    vinit = ( ones(size(x))+0.3*noise ) * p.v0;
    % vinit = u0 + v0 - uinit; 
    
    init = [uinit; vinit];
    
end