% Programa principal para resolver la ecuación: 
% 
%       (k(x)u')'  = f(x) en el intervalo (a,b)
% Con condiciones de Dirchlet o Robin dadas por:
%       alpha0*u(0) - beta0*u'(0) = C0;
%       alpha1*u(1) + beta1*u'(1) = C1;
%  
% con  un esquema mimético.
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
%   . La solución exacta esta dada por: 
%         u(x) = -(11/6)x^3 + (15/2)x^2 - (20/3)x + 1; 
%  El código esta pensando para introducir los Miméticos en un curso corto
%  usando la libreria MOLE  
%                GECS   11/10/2024
%-------------------------------------------------------------
clc, clear
close all

%addpath('../mole_MATLAB')
addpath('C:\Users\TATUY\Dropbox\Mimetico2024\MOLE\mole_MATLAB')
west = 0;  % Domain's limits
east = 2;
C0 = 1; % West BC 
C1 = 3; % East BC
fuente = @(x)-33/2*x.^2 + 30*x - 20/3; %f(x)
SOL    = @(x)-(11/6)*x.^3 + (15/2)*x.^2 - (20/3)*x + 1; % Sol exacta
K      = @(x)x;
nelem = 40; % number of cells
dx =(east-west)/nelem;% Step size
orden = 2; % Order of accuracy
xgridSca = [west west+.5*dx: dx :east-.5*dx east]; % 1D Staggered grid
d =  west+.5*dx: dx :east-.5*dx;  % Staggered grid sin frontera
xgridVec = west : dx : east;
G = grad(orden,nelem,dx);
D = div(orden,nelem,dx);
BC = robinBC_GECS(orden,nelem,dx,1,0);
f = fuente(d);
F = [C0 f C1];
Z= diag(K(xgridVec),0);
KK = BC + D*(Z*G);
UDM = (KK\F'); % Calculo de la Solución Usando Mimeticos.


% Plot result
plot(xgridSca, SOL(xgridSca))
hold on
plot(xgridSca, UDM, 'o')
legend('Sol Exacta', 'Sol aproximada', 'Location', 'NorthEast')
title('Ejemplo 3: Solución con BC Dirichlet')
xlabel('x')
ylabel('u(x)')


