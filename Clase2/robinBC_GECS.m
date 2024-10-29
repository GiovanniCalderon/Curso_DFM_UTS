function BC = robinBC_GECS(k, m, dx, a,b)
% Returns a m+2 by m+2 one-dimensional mimetic boundary operator that 
% imposes a boundary condition of Robin's type
%
% Parameters:
%                k : Order of accuracy
%                m : Number of cells
%               dx : Step size
%                a : Dirichlet Coefficient
%                b : Neumann Coefficient

    A = sparse(m+2, m+2);
    A(1, 1) = a;
    A(end, end) = a;
    
    B = sparse(m+2, m+1);
    B(1, 1) = -b;
    B(end, end) = b;
    
    G = grad(k, m, dx);
    %%%%%%%%%%%%% Bloque añadido
    if b~=0
       % B=mimeticB(k,m);

        %B=0*B; 

        B = sparse(m+2, m+1);
        B(1,1)=-1;B(2,1)=1/8;B(2,2)=-1/8;
        B(3,1)=-1/8;B(3,2)=1/8;

        B(end,end)=1; B(end-1,end-1)=1/8; B(end-1,end)=-1/8;
        B(end-2,end-1)=-1/8; B(end-2,end)=1/8;
    end
    %%%%%%%%%%%%% 
    BC = A + b*B*G;
end
