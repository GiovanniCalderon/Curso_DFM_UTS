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
%addpath('C:\Users\UIS\Dropbox\Mimetico2024\MOLE\mole_MATLAB')
addpath('C:\Users\TATUY\Dropbox\Mimetico2024\MOLE\mole_MATLAB')
west = 0; % Domain's limits
east = 1; 
M = [10 40 80 120 160 200 400 600 1200]; 
dim = length(M);
k = 1 ; v = 1;
fuente = @(x)v*(10*sin(pi*x.^2) - (pi*sin((pi*x)/2))/2 + 20*x.^2*pi.*cos(pi*x.^2)) + ...
         (k*pi*(pi*cos((pi*x)./2) - 240*x.*cos(pi*x.^2) + 160*x.^3.*pi.*sin(pi*x.^2)))./4;

SOL =  @(x)cos(.5*pi*x) + 10*x.*sin(pi*x.^2);
TipoFrontera   = 'Robin';% 'Dirichlet'; % 
xx = west:.0001:east;    %solo para graficar
H = zeros(dim,1);
ErrorNorMaxMM = zeros(dim,1);
ErrorNorL2MM  = zeros(dim,1);
Pe = zeros(dim,1);
for pasos = 1:dim
    nelem = M(pasos); % number of cells
    fprintf('----- Calculando para %4i elemetos \n',nelem);
    dx =(east-west)/nelem;% Step size
    Pe(pasos) = abs(v)*dx/(2*k); % Peclet % en algunos dividen por 2 tambien.
    nnode = nelem+1;
    orden = 2; % Order of accuracy
    L = lap(orden,nelem,dx); % 1D Mimetic laplacian operator
    G = grad(orden,nelem,dx);
    D = div(orden,nelem,dx);
    ID = interpol(nelem, .5);
    xgridSca = [west west+.5*dx: dx :east-.5*dx east]; % 1D Staggered grid
    d =  west+.5*dx: dx :east-.5*dx ;
    %d = linspace(west+.5*dx,east-.5*dx,nelem); % Staggered grid sin frontera
    xgridVec = linspace(west,east,nnode);
    UDM1 = zeros(nnode+1,1); % matriz solucion: METODO MIMETICO
    % =========================================================================
    f = fuente(d);       
    if strcmp(TipoFrontera,'Dirichlet')
        alpha1 = 1;   alpha2 = 1; 
        beta1  = 0; beta2  = 0;
        C0 = 1;
        C1 = 0;

      % Mimetic Method  
        BC = robinBC_GECS(orden,nelem,dx,alpha1,0);
        %BC = robinBC_GECSv2(orden,nelem,dx,alpha1,alpha2,beta1,beta2);
        F = [C0 f C1];
        KK = BC - k*D*G + v*D*ID;
        UDM1 = (KK\F'); 
      % End Mimetic method.
    end

    if strcmp(TipoFrontera,'Robin')
        alpha1 = 1; alpha2 = 1;
        beta1  = 1; beta2  = 1;
        C0 = 1;
        C1 = -20.5*pi;
        F = [C0 f C1];
      % Calculo de la Solución Usando Mimeticos.
        BC = robinBC_GECSv2(orden,nelem,dx,alpha1,alpha2,beta1,beta2);
        %BC = robinBC_GECS(orden,nelem,dx,alpha1,1);
        KK = BC - k*D*G + v*D*ID;
        UDM1 = (KK\F'); % Calculo de la Solución Usando Mimeticos.
    end
    % POSTPROCESO DEL METODO:
    Uex_DM = SOL(xgridSca); % Eval sol exacta con DM:
    % NORMA L1
    ErrorNorMaxMM(pasos) = dx*sum(abs(Uex_DM' - UDM1));   
    H(pasos) = dx;
  
    HH = H(pasos)*ones(nelem+2,1);  
    HH(1) = .5*HH(1);   HH(end) = .5*HH(end);
    %ErrorNorL2MM(pasos)  = sqrt(sum((HH(pasos).*(Uex_DM'-UDM1)).^2));
    ErrorNorL2MM(pasos)  = sqrt(sum(HH(pasos).*(Uex_DM'-UDM1).^2));
    ErrorNorMaxMM(pasos) = sum(HH(pasos).*abs(Uex_DM' - UDM1));
    if strcmp(TipoFrontera,'Dirichlet')
       figure(1)
       plot(xx,SOL(xx),'-k',xgridSca,UDM1,'r')
       legend('Exacta','Mimetico')

       figure(2)
       plot(xgridSca,abs(Uex_DM' - UDM1),'r')
       legend('Mimetico')
       title('Error')
    end 
    if strcmp(TipoFrontera,'Robin')
       figure(1)
       plot(xx,SOL(xx),'-k',xgridSca,UDM1,'r')
       legend('Exacta','Mimetico')

       figure(2)
       plot(xgridSca,abs(Uex_DM' - UDM1),'r')
       legend('Mimetico')
       title('Magnitud del Error puntual')
    end

    pause   
end
    
% POSTPROCESO DE DATOS:
fprintf('----- \b\b\b\033[1mRESULTADOS\033[0m -----  \n');
for pasos = 1:dim
    nelem = M(pasos);
    fprintf('----- Numero de elemetos = %4i \n',nelem);
    fprintf('----- Numero de Peclet   = %6.2e \n',Pe(pasos));
    fprintf('----- Error en norma l_1:  %8.4e\n',ErrorNorMaxMM(pasos))
    fprintf('----- Error en norma l_2:  %8.4e\n',ErrorNorL2MM(pasos))
    if pasos > 1
        p1 = log(ErrorNorMaxMM(pasos)/ErrorNorMaxMM(pasos-1))/log(H(pasos)/H(pasos-1));
        p2 = log(ErrorNorL2MM(pasos)/ErrorNorL2MM(pasos-1))/log(H(pasos)/H(pasos-1));
        fprintf('----- Orden C.  en norma l_1:  %10.4e\n',p1)
        fprintf('----- Orden C.  en norma l_2:  %10.4e\n',p2)
    end  
    fprintf(' \n')
end
fprintf('----------------------  \n');


if strcmp(TipoFrontera,'Dirichlet')  
    figure(10)
    loglog(H,ErrorNorMaxMM,'-dr',H,ErrorNorL2MM,'-db')
    legend('Orden en l_1','Orden en l_2')
    title('Orden de Convergencia')
    
end
if strcmp(TipoFrontera,'Robin')
    figure(20)
    loglog(H,ErrorNorMaxMM,'-dr',H,ErrorNorL2MM,'-db')
    legend('Orden en l_1','Orden en l_2')
    title('Orden de Convergencia')
end


  
  

