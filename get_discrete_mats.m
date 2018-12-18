function [F, G, H, omega] = get_discrete_mats(x, u, L, delta_T)
    A = [0 0 -u(1)*sin(x(3)) 0 0 0; 
         0 0 u(1)*cos(x(3)) 0 0 0; 
         0 0 0 0 0 0; 
         0 0 0 0 0 -u(3)*sin(x(6)); 
         0 0 0 0 0 u(3)*cos(x(6)); 
         0 0 0 0 0 0];
    B = [cos(x(3)) 0 0 0; 
         sin(x(3)) 0 0 0; 
         tan(u(2))/L u(1)/(L*(cos(u(2)))^2) 0 0;
         0 0 cos(x(6)) 0; 
         0 0 sin(x(6)) 0; 
         0 0 0 1];
    
    x52 = x(5)-x(2); x41 = x(4)-x(1); x14 = -x41; x25 = -x52;
    C = [x52/(x41^2+x52^2) x14/(x41^2+x52^2) -1 x25/(x41^2+x52^2) x41/(x41^2+x52^2) 0; 
         x14/sqrt(x14^2+x25^2) x25/sqrt(x14^2+x25^2) 0 x41/sqrt(x14^2+x25^2) x52/sqrt(x14^2+x25^2) 0; 
         x52/(x14^2+x25^2) x14/(x14^2+x25^2) 0 x25/(x14^2+x25^2) x41/(x14^2+x25^2) -1; 
         0 0 0 1 0 0; 
         0 0 0 0 1 0];
    
    F = eye(size(A)) + delta_T.*A;
    G = delta_T.*B;
    H = C;
    omega = delta_T.*eye(size(A));
end