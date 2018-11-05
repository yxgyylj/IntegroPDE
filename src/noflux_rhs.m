%% define the discretized ODE system for periodic BC

function out = noflux_rhs(t,n,p)

    if t > 30 && p.flag == 1
        p.k2 = .432;
    end
    
    p.Nx = length(p.x); 
    u = n(1:p.Nx); 
    v = n(p.Nx+1:end); 

    int_uu = zeros(p.Nx,1);
    int_uv = zeros(p.Nx,1);
    int_vu = zeros(p.Nx,1);
    
    % Nonlocal Integral for u and v equations
    uM = (p.u0+p.v0); 
    for i=1:p.Nx
        intgrd_uu = Gr(p.x-p.x(i),p.r1).*u;%.*(1-u/uM);
        intgrd_uv = Gr(p.x-p.x(i),p.r2).*v;
        intgrd_vu = Gr(p.x-p.x(i),p.r2).*u;
        int_uu(i) = sum(intgrd_uu)*p.dx; 
        int_uv(i) = sum(intgrd_uv)*p.dx; 
        int_vu(i) = sum(intgrd_vu)*p.dx;
    end                                                                 
    
    % flux terms for u and v equations at half grids
    fluxu = (u(2:end) - u(1:end-1))*p.D1/p.dx ... % linear term
        + (u(1:end-1) + u(2:end)).*(u(2:end) - u(1:end-1) ...
        + v(2:end) - v(1:end-1))*p.eta1/(2*p.dx) ... % nonlinear term
        - (u(1:end-1).*(int_uu(1:end-1)*p.k1 + int_uv(1:end-1)*p.k2)...
        +u(2:end).*(int_uu(2:end)*p.k1 + int_uv(2:end)*p.k2))/2/p.mu1; % nonlocal term
    fluxv = (v(2:end) - v(1:end-1))*p.D2/(2*p.dx) ... % linear term
        + (v(1:end-1) + v(2:end)).*(u(2:end) - u(1:end-1) ...
        + v(2:end) - v(1:end-1))*p.eta2/(2*p.dx) ... % nonlinear term
        - (v(1:end-1).*int_vu(1:end-1)...
        + v(2:end).*int_vu(2:end))*p.k2/2/p.mu2; % nonlocal term
    
    % shift vectors for flux terms, which should satisfy B.Cs and 
    % conservation of mass
    fluxum = [-fluxu(1); fluxu];
    fluxup = [fluxu; -fluxu(end)];
    fluxvm = [-fluxv(1); fluxv];
    fluxvp = [fluxv; -fluxv(end)];
        
    % Central Difference Approximation for u and v equations
    du = (fluxup - fluxum)/(p.dx);
    du = du(:);
    dv = (fluxvp - fluxvm)/(p.dx); 
    dv = dv(:); 
    
    out = [du; dv]; 
    
end
