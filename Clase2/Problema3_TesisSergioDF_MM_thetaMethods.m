% Programa principal para resolver la ecuación de 
% 
%       -ku'' +vu' = f(x) en el intervalo (0,1)
% Con condiciones de Dirchlet o Robin dadas por:
%       alpha0*u(0) - beta0*u'(0) = C0;
%       alpha1*u(1) + beta1*u'(1) = C1;
%   k, v : constantes v>0
% con esquemas en DF: Upwind 1er, DC de 2do y SG, mas un esquema mimético.
% PARAMETROS NECESARIOS :
%    nelem : Número de elementos de la malla
%    [a,b] : Valores extremos del dominio de definicion
%    fuente: Nombre de la funcion que define el termino fuente del problema
%    TipoFrontera: Se tienen definidas dos tipos de condiciones de frontera
%                  'Dirichlet' y 'Robin'. 
% Observación: 
%   . Se usa la libreria MOLE de Jose Castillo
%   . Se lleva un esquema de DF de 2do orden y 2do orden centrado para el 
%     termino convectivo para comparar la convergencia.
%   . Para este ejemplo vamos a tomar k = , v=  distinto de cero. 
%   .  La solución exacta esta dada por: 
%         u(x) = cos(.5*pi*x)+10*x.*sin(pi*x.^2); 
%  El código esta pensando para el problema modelo de la tesis de Sergio  
%                GECS   12/08/2024
%-------------------------------------------------------------
clc, clear, close all
format long
addpath('C:\Users\UIS\Dropbox\Mimetico2024\MOLE\mole_MATLAB')
addpath('C:\Users\TATUY\Dropbox\Mimetico2024\MOLE\mole_MATLAB')
west = 0; % Domain's limits
east = 1;
%M = [10 20 40]; 
M = [40 80 120 160 200 400 600]; 
%M = [500 700 1000 1500]; 
dim = length(M);
k = 150 ; v = 1;
fuente = @(x)v*(10*sin(pi*x.^2) - (pi*sin((pi*x)/2))/2 + 20*x.^2*pi.*cos(pi*x.^2)) + ...
         (k*pi*(pi*cos((pi*x)./2) - 240*x.*cos(pi*x.^2) + 160*x.^3.*pi.*sin(pi*x.^2)))./4;


SOL =  @(x)cos(.5*pi*x)+10*x.*sin(pi*x.^2);
TipoFrontera   = 'Robin';% 'Dirichlet'; % 
xx = west:.0001:east;    %solo para graficar
H = zeros(dim,1);
ErrorNorMaxMM  = zeros(dim,1);
ErrorNorMaxDF1 = zeros(dim,1);
ErrorNorMaxDF2 = zeros(dim,1);
ErrorNorMaxDF3 = zeros(dim,1);
ErrorNorMaxDF4 = zeros(dim,1);

ErrorNorL2MM  = zeros(dim,1);
ErrorNorL2DF1 = zeros(dim,1);
ErrorNorL2DF2 = zeros(dim,1);
ErrorNorL2DF3 = zeros(dim,1);
ErrorNorL2DF4 = zeros(dim,1);
for pasos = 1:dim
    nelem = M(pasos); % number of cells
    fprintf('----- RESULTADOS -----  \n');
    fprintf('----- Numero de elemetos = %4i \n',nelem);
    dx =(east-west)/nelem;% Step size
    Pe = abs(v)*dx/(2*k); % Peclet % en algunos dividen por 2 tambien.
    fprintf('----- Numero de Peclet = %10.2e \n',Pe);
    nnode = nelem+1;
    orden = 2; % Order of accuracy
    L = lap(orden,nelem,dx); % 1D Mimetic laplacian operator
    G = grad(orden,nelem,dx);
    D = div(orden,nelem,dx);
    II = interpolD(nelem, .5);
    ID = interpol(nelem, .5);
    xgridSca = [west west+.5*dx: dx :east-.5*dx east]; % 1D Staggered grid
    d =  west+.5*dx: dx :east-.5*dx ;
    %d = linspace(west+.5*dx,east-.5*dx,nelem); % Staggered grid sin frontera
    xgridVec = linspace(west,east,nnode);
    UDF1 = zeros(nnode,1);   % matriz solucion: DIFERENCIAS FINITAS Theta Methods
    UDF2 = zeros(nnode,1);   % matriz solucion: DIFERENCIAS FINITAS Dif Centradas orden 2
    UDF3 = zeros(nnode,1);   % matriz solucion: DIFERENCIAS FINITAS UPWING_O2
    UDF4 = zeros(nnode,1);   % matriz solucion: DIFERENCIAS FINITAS QUICK
    UDM1 = zeros(nnode+1,1); % matriz solucion: METODO MIMETICO
    UDM2 = zeros(nnode+1,1); % matriz solucion: METODO MIMETICO
    % =========================================================================
    I = diag(ones(nelem+2,1));
    I1= diag(zeros(nelem+2,1)); I1(1,1)=1;I1(end,end)=1;
    f = fuente(d);       
    if strcmp(TipoFrontera,'Dirichlet')
        alpha0 = 1;
        %alpha1 = 1;
        C0 = 1;
        C1 = 0;
        aa = ones(nnode-2,1);
      % Inicio DF Centradas
        W = sparse(diag(-(1+Pe)*aa(2:nnode-2),-1)) + sparse(2*diag(aa,0))+sparse((Pe-1)*diag(aa(1:nnode-3),1));
        S = (dx^2/k)*fuente(xgridVec(2:end-1));
        S(1) = S(1) +(1+Pe)*C0;
        S(nnode-2) = S(nnode-2) + (1-Pe)*C1;
        UDF2(2:end-1) = W\S';        
        UDF2(1) = C0;
        UDF2(end) = C1;
      % end DF Centradas

      % Theta Methods DF
        s = k/dx^2; c = v/dx;
       % phi = v*dx/(3*k);    
        phi = 2*(k-2)*dx/(v*dx^2+12);
        a1 = -(s+phi*c+.5*c*(1-phi)); a2 = 2*s+c*phi; a3 = (-s+.5*c*(1-phi));
        W = sparse(diag(a1*aa(2:nnode-2),-1)) + sparse(a2*diag(aa,0))+sparse(a3*diag(aa(1:nnode-3),1));
        S = fuente(xgridVec(2:end-1));
        S(1) = S(1) - a1*C0;
        S(nnode-2) = S(nnode-2) - a3*C1;
        UDF1(2:end-1) = W\S';        
        UDF1(1) = C0;
        UDF1(end) = C1;
      % end Theta Methods DF

      % Upwind_O2 Methods DF
        s = k/dx^2; c = v/dx;
        a1 = .5*c; a2 = -(s+2*c); a3 = 2*s+1.5*c; a4 = -s;
        W = sparse(diag(a1*aa(3:nnode-2),-2)) + sparse(diag(a2*aa(2:nnode-2),-1)) +...
            sparse(a3*diag(aa,0)) + sparse(a4*diag(aa(1:nnode-3),1));
        W(1,1) = 2*s; W(1,2) = .5*c-s;
        S = fuente(xgridVec(2:end-1));
        S(1) = S(1) + (s+.5*c)*C0;
        S(2) = S(2) -.5*c*C0;
        S(nnode-2) = S(nnode-2) + s*C1;
        UDF3(2:end-1) = W\S';        
        UDF3(1) = C0;
        UDF3(end) = C1;
      % end Upwind_O2 Methods DF


      % Quick Methods DF
        s = k/dx^2; c = v/dx;
        a1 = (1/8)*c; a2 = -(s+(7/8)*c); a3 = 2*s+(3/8)*c; a4 = (3/8)*c-s;
        W = sparse(diag(a1*aa(3:nnode-2),-2)) + sparse(diag(a2*aa(2:nnode-2),-1)) +...
            sparse(a3*diag(aa,0)) + sparse(a4*diag(aa(1:nnode-3),1));
        W(1,1) = 2*s-(3/8)*c; W(1,2) = -(s-(7/8)*c); W(1,3) = -(1/8)*c;
        S = fuente(xgridVec(2:end-1));
        S(1) = S(1) + (s+(3/8)*c)*C0;
        S(2) = S(2) -(1/8)*c*C0;
        S(nnode-2) = S(nnode-2) -((3/8)*c-s)*C1;
        UDF4(2:end-1) = W\S';        
        UDF4(1) = C0;
        UDF4(end) = C1;
      % end Quick Methods DF

      % Theta Methods MM     
        phiD = 2*(k-2)*dx/(v*dx^2+12);%.5*(1-c);
        DD = zeros(nnode+1,1);    
        DD(2) = 8*(phiD-1)/(3*dx)*UDM2(1)+(1/dx)*(3-4*phiD)*UDM2(2)+1/(3*dx)*(4*phiD-1)*UDM2(3);
        DD(nnode) = 1/(3*dx)*(4*phiD-3)*UDM2(nnode-1) + (1/dx)*(1-4*phiD)*UDM2(nnode) + 1/(3*dx)*(8*phiD)*UDM2(nnode+1);
        for j = 3:nnode-1
            DD(j) = (phiD/dx)*(UDM2(j+1)-UDM2(j))+((1-phiD)/dx)*(UDM2(j)-UDM2(j-1));
        end
        BC = robinBC_GECS(orden,nelem,dx,alpha0,0);
        F = [C0 f C1];
        KK = BC - k*D*G + v*DD;
        UDM2 = (KK\F');
      % end Theta Methods MM

      % Mimetic Method    
        BC = robinBC_GECS(orden,nelem,dx,alpha0,0);
        F = [C0 f C1];
         KK = BC - k*D*G + v*D*ID;
        UDM1 = (KK\F'); 
      % End Mimetic method.

    end



    if strcmp(TipoFrontera,'Robin')
        alpha0 = 1; alpha1 = 1;
        beta0  = 1; beta1  = 1;
        C0 = 1;
        C1 = -20.5*pi;
        % UPWIND Orden 1
        % a = -(1+2*Pe)*ones(nnode,1); b = (2+2*Pe)*ones(nnode,1); c = ones(nnode,1);
        % W = sparse(diag(a(2:nnode-2),-1)) + sparse(diag(b,0))-sparse(diag(c(1:nnode-3),1));
        % W(1,1) = (1+2*Pe)*alpha0+2+2*Pe;
        % W(1,2) = 2+2*Pe;
        % W(nelem+1,nelem+1)= 2+2*Pe+2*dx*alpha1/beta1;
        % W(nelem+1,nelem) = -(2+2*Pe);
        % S = (dx^2/k)*fuente(xgridVec);
        % S(1) = S(1) +(1+2*Pe)*2*dx*C0/beta0;
        % S(nelem+1) = S(nelem+1) + 2*dx*C1/beta1;
        % UDF1 = W\S';      
      % end de UPWIND Orden 1

      % Inicio DF Centradas
        oo = ones(nelem+1,1);
        a = -k/dx^2-v/(2*dx);   b= 2*k/dx^2;   c= v/(2*dx)-k/dx^2;
        aa = a*oo; bb = b*oo; cc = c*oo;
        W = sparse(diag(aa(2:nelem+1),-1)) +1*sparse(diag(bb,0)) + sparse(diag(cc(1:nelem),1));
        W(1,1) = -a*2*dx*alpha0/beta0 + b;
        W(1,2) = a+c;
        W(nelem+1,nelem+1)= b-c*alpha1*2*dx/beta1;
        W(nelem+1,nelem) = a+c;
        S = fuente(xgridVec);
        S(1) = S(1) - a*2*dx*C0/beta0;
        S(nelem+1) = S(nelem+1) - c*2*dx*C1/beta1;
        UDF2 = W\S';
        % end DF Centradas

      % Theta Methods DF
        s = k/dx^2; c = v/dx;
       % phi = v*dx/(3*k);    
        phi = 2*(k-2)*dx/(v*dx^2+12);
        a1 = -(s+phi*c+.5*c*(1-phi)); a2 = 2*s+c*phi; a3 = (-s+.5*c*(1-phi));
        W = sparse(diag(a1*oo(2:nelem+1),-1)) + sparse(a2*diag(oo,0))+sparse(a3*diag(oo(1:nelem),1));
        W(1,1) = -a1*2*dx*alpha0/beta0 + a2;
        W(1,2) = a1+a3;
        W(nelem+1,nelem+1)= a2-a3*alpha1*2*dx/beta1;
        W(nelem+1,nelem) = a1+a3;
        S = fuente(xgridVec);
        S(1) = S(1) - a1*2*dx*C0/beta0;
        S(nelem+1) = S(nelem+1) - a3*2*dx*C1/beta1;
        UDF1 = W\S';
      % end Theta Methods DF
      
      % Upwind O(2) Methods DF
        s = k/dx^2; c = v/dx;
        a1 = .5*c; a2 = -(s+2*c); a3 = 2*s+1.5*c; a4 = -s;
        W = sparse(diag(a1*oo(3:nnode),-2)) + sparse(diag(a2*oo(2:nnode),-1)) +...
            sparse(a3*diag(oo,0)) + sparse(a4*diag(oo(1:nnode-1),1));
        W(1,1) = alpha0+1.5*beta0/dx; W(1,2) = -2*beta0/dx; W(1,3) = beta0/(2*dx);
        W(2,1) = -(s+.5*c); W(2,2) = 2*s; W(2,3) = .5*c-s;
        W(nnode,nnode-2) = beta1/(2*dx); W(nnode,nnode-1)= -2*beta1/dx; W(nnode,nnode) = alpha1+1.5*beta1/dx;
        S = fuente(xgridVec);
        S(1) = C0;
        S(nnode) = C1;
        UDF3 = W\S';        
      % end Upwind O(2) Methods DF

      % Quick Methods DF
        s = k/dx^2; c = v/dx;
        a1 = (1/8)*c; a2 = -(s+(7/8)*c); a3 = 2*s+(3/8)*c; a4 = (3/8)*c-s;
        W = sparse(diag(a1*oo(3:nnode),-2)) + sparse(diag(a2*oo(2:nnode),-1)) +...
            sparse(a3*diag(oo,0)) + sparse(a4*diag(oo(1:nnode-1),1));
        W(1,1) = alpha0+1.5*beta0/dx; W(1,2) = -2*beta0/dx; W(1,3) = beta0/(2*dx);
        W(2,1) = -(s+(3/8)*c); W(2,2) = 2*s-(3/8)*c; W(2,3) = -s+(7/8)*c; W(2,4) = -(1/8)*c;
       % W(2,1) = -(s+(1/2)*c); W(2,2) = 2*s; W(2,3) = -s+.5*c; W(2,4) = 0;
        W(nnode,nnode-2) = beta1/(2*dx); W(nnode,nnode-1)= -2*beta1/dx; W(nnode,nnode) = alpha1+1.5*beta1/dx;
        S = fuente(xgridVec);
        S(1) = C0;
        S(nnode) = C1;
        UDF4 = W\S';        
      % end Quick Methods DF

      % Calculo de la Solución Usando Mimeticos.
        BC = robinBC_GECS(orden,nelem,dx,alpha0,1);
        %F = fuente(xgridSca); F = [C0+F(1) f C1+F(end)];
        F = [C0 f C1];
        %KK = BC - k*D*G + v*ID*G;
         KK = BC - k*D*G + v*D*ID;
        UDM1 = (KK\F'); % Calculo de la Solución Usando Mimeticos.
    end
  % POSTPROCESO DEL METODO:
    Uex_DM = SOL(xgridSca); % Eval sol exacta con DM:
    Uex_DF =  SOL(xgridVec); % Eval sol exacta con DF:
  % NORMA max: ||e||_inf
    %ErrorNorMaxMM(pasos) = max(abs(Uex_DM' - UDM1));
    %ErrorNorMaxDF(pasos)  = max(abs(Uex_DF' - UDF1));
% NORMA L1
    ErrorNorMaxMM(pasos)   = dx*sum(abs(Uex_DM' - UDM1));   
    ErrorNorMaxDF1(pasos)  = dx*sum(abs(Uex_DF' - UDF1));
    ErrorNorMaxDF2(pasos)  = dx*sum(abs(Uex_DF' - UDF2));
    ErrorNorMaxDF3(pasos)  = dx*sum(abs(Uex_DF' - UDF3));
    ErrorNorMaxDF4(pasos)  = dx*sum(abs(Uex_DF' - UDF4));
    H(pasos) = dx;
  
%    ErrorNorL2MM(pasos) = sqrt(H(pasos)*sum((Uex_DM'-UDM1).^2));
%    ErrorNorL2DF(pasos) = sqrt(H(pasos)*sum((Uex_DF'-UDF1).^2));

    HH = H(pasos)*ones(nelem+2,1);  
    HH(1) = .5*HH(1);   HH(end) = .5*HH(end);
    ErrorNorL2MM(pasos)  = sqrt(sum((HH(pasos).*(Uex_DM'-UDM1)).^2));
    ErrorNorL2DF1(pasos) = sqrt(sum((H(pasos)*(Uex_DF'-UDF1)).^2));
    ErrorNorL2DF2(pasos) = sqrt(sum((H(pasos)*(Uex_DF'-UDF2)).^2));
    ErrorNorL2DF3(pasos) = sqrt(sum((H(pasos)*(Uex_DF'-UDF3)).^2));
    ErrorNorL2DF4(pasos) = sqrt(sum((H(pasos)*(Uex_DF'-UDF4)).^2));
% ERROR LOCAL normalizado
    kk1 = abs(Uex_DM'-UDM1)./abs(Uex_DM');
    kk2 = abs(Uex_DF'-UDF1)./abs(Uex_DF');

    %ErrorNorL2MM(pasos) = (1/(nnode+1))*sum(kk1(2:end-1));
    %ErrorNorL2DF(pasos) = (1/nnode)*sum(kk2(2:end-1));

    % ErrorNorL2MM(pasos)  = sqrt(sum((H(pasos)*(Uex_DM'-UDM1)).^2));
    % ErrorNorL2DF1(pasos) = sqrt(sum((H(pasos)*(Uex_DF'-UDF1)).^2));
    % ErrorNorL2DF2(pasos) = sqrt(sum((H(pasos)*(Uex_DF'-UDF2)).^2));
    % ErrorNorL2DF3(pasos) = sqrt(sum((H(pasos)*(Uex_DF'-UDF3)).^2));
    % ErrorNorL2DF4(pasos) = sqrt(sum((H(pasos)*(Uex_DF'-UDF4)).^2));
    
    if strcmp(TipoFrontera,'Dirichlet')
       figure(1)
       plot(xx,SOL(xx),'-k',xgridSca,UDM1,'r',xgridVec,UDF1,'b',xgridVec,UDF2,'g',xgridVec,UDF3,'m',xgridVec,UDF4,'c')
       legend('Exacta','Mimetico', 'DF Theta', 'DF centradas','DF Upwing-O2','DF Quick')

       figure(2)
       plot(xgridSca,abs(Uex_DM' - UDM1),'r',xgridVec,abs(Uex_DF' - UDF1),'b',...
            xgridVec,abs(Uex_DF' - UDF2),'g',xgridVec,abs(Uex_DF' - UDF3),'m',xgridVec,abs(Uex_DF' - UDF4),'c')
       legend('Mimetico','DF Theta', 'DF centradas','DF Upwing-O2','DF Quick')
    end 
    if strcmp(TipoFrontera,'Robin')
       figure(1)
       plot(xx,SOL(xx),'-k',xgridSca,UDM1,'r',xgridVec,UDF1,'b',xgridVec,UDF2,'g',xgridVec,UDF3,'m',xgridVec,UDF4,'c')
       legend('Exacta','Mimetico', 'DF Theta', 'DF centradas','DF Upwing-O2','DF Quick')

       figure(2)
       plot(xgridSca,abs(Uex_DM' - UDM1),'r',xgridVec,abs(Uex_DF' - UDF1),'b',...
            xgridVec,abs(Uex_DF' - UDF2),'g',xgridVec,abs(Uex_DF' - UDF3),'m',xgridVec,abs(Uex_DF' - UDF4),'c')
       legend('Mimetico','DF Theta', 'DF centradas','DF Upwing-O2','DF Quick')
    end

    fprintf('-----  Errores: norma l_1               norma l_2 \n')
    fprintf('-----  MDM   :  %14.4e       %14.4e\n',ErrorNorMaxMM(pasos),ErrorNorL2MM(pasos))
    fprintf('-----  T-DF  :  %14.4e       %14.4e\n',ErrorNorMaxDF1(pasos),ErrorNorL2DF1(pasos))
    fprintf('-----  DFC   :  %14.4e       %14.4e\n',ErrorNorMaxDF2(pasos),ErrorNorL2DF2(pasos))
    fprintf('-----  Upwind:  %14.4e       %14.4e\n',ErrorNorMaxDF3(pasos),ErrorNorL2DF3(pasos))
    fprintf('-----  Quick :  %14.4e       %14.4e\n',ErrorNorMaxDF4(pasos),ErrorNorL2DF4(pasos))
    fprintf('----------------------  \n');
    pause   

end


if strcmp(TipoFrontera,'Dirichlet')
    figure(10)
    loglog(H,ErrorNorMaxMM,'-dr',H,ErrorNorMaxDF1,'-sb',H,ErrorNorMaxDF2,'-og',H,ErrorNorMaxDF3,'-vm',H,ErrorNorMaxDF4,'->c')
    legend('Mimetico','DF Theta', 'DF centradas','DF Upwing-O2','DF Quick')
    title('Error en norma l_1')

    figure(20)
    loglog(H,ErrorNorL2MM,'-dr',H,ErrorNorL2DF1,'-sb',H,ErrorNorL2DF2,'-og',H,ErrorNorL2DF3,'-vm',H,ErrorNorL2DF4,'->c')
    legend('Mimetico','DF Theta', 'DF centradas','DF Upwing-O2','DF Quick')
    title('Error en norma l_2')
end
if strcmp(TipoFrontera,'Robin')
    figure(10)
    loglog(H,ErrorNorMaxMM,'-dr',H,ErrorNorMaxDF1,'-sb',H,ErrorNorMaxDF2,'-og',H,ErrorNorMaxDF3,'-vm',H,ErrorNorMaxDF4,'->c')
    legend('Mimetico','DF Theta', 'DF centradas','DF Upwing-O2','DF Quick')
    title('Error en norma l_1')

    figure(20)
    loglog(H,ErrorNorL2MM,'-dr',H,ErrorNorL2DF1,'-sb',H,ErrorNorL2DF2,'-og',H,ErrorNorL2DF3,'-vm',H,ErrorNorL2DF4,'->c')
    legend('Mimetico','DF Theta', 'DF centradas','DF Upwing-O2','DF Quick')
    title('Error en norma l_2')
end


  
  

