%%------------- 4-th order Runge–Kutta method for ODEs --------------------
%   input:
%       rhs_fun :   right hand side function, depents on time t and a vector u 
%       t       :   current time
%       u       :   input vector, the function you want to solve
%       dt      :   increatment in time, default to be 0.01
%
%   written by: Xige Yang
%%-------------------------------------------------------------------------

function out = RK4(rhs_fun,t,u,dt)
    if nargin == 3
        dt = .01;
    end
    
    k1 = dt*rhs_fun(t,u);
    tmp = u + 0.5*k1;
    k2 = dt*rhs_fun(t+0.5*dt,tmp);
    tmp = u + 0.5*k2;
    k3 = dt*rhs_fun(t+0.5*dt,tmp);
    tmp = u + k3;
    k4 = dt*rhs_fun(t+dt,tmp);
    out = u+1/6*(k1+k4)+1/3*(k2+k3);
end


