clc
close all
clear variables

grado_polinomio = 12;
nodos_int = grado_polinomio +1;

[chis, ws] = Integracion_GLL(nodos_int+1);
[chis_NE, ws_NE] = Integracion_GLL(nodos_int);

PHI_NL = Base_NodalLagrange(grado_polinomio);
dPHI_NL = zeros(size(PHI_NL));
for n = 1:nodos_int
    dPHI_NL(n, 2:end) = polyder(PHI_NL(n, :));
end
PP_NL = zeros(nodos_int);
PdP_NL = size(PP_NL);
dPdP_NL = size(PP_NL);

PHI_NE = Base_NodalEspectral(grado_polinomio);
dPHI_NE = zeros(size(PHI_NE));
for n = 1:nodos_int
    dPHI_NE(n, 2:end) = polyder(PHI_NE(n, :));
end
PP_NE = size(PP_NL);
PdP_NE = size(PP_NL);
dPdP_NE = size(PP_NL);

PHI_M = Base_Modal(grado_polinomio, 1, 1);
dPHI_M = zeros(size(PHI_M));
for n = 1:nodos_int
    if n<nodos_int
        dPHI_M(n, end-n+1:end) = polyder(PHI_M(n, :));
    else
        dPHI_M(n, end) = polyder(PHI_M(n, :));
    end
end
PP_M = size(PP_NL);
PdP_M = size(PP_NL);
dPdP_M = size(PP_NL);

for i = 1:nodos_int
    for j = 1:nodos_int
        PP_NL(i, j) = sum(ws.*polyval(conv(PHI_NL(i, :), PHI_NL(j, :)), chis));
        PdP_NL(i, j) = sum(ws.*polyval(conv(PHI_NL(i, :), dPHI_NL(j, :)), chis));
        dPdP_NL(i, j) = sum(ws.*polyval(conv(dPHI_NL(i, :), dPHI_NL(j, :)), chis));
        
        PP_NE(i, j) = sum(ws_NE.*polyval(conv(PHI_NE(i, :), PHI_NE(j, :)), chis_NE));
        PdP_NE(i, j) = sum(ws_NE.*polyval(conv(PHI_NE(i, :), dPHI_NE(j, :)), chis_NE));
        dPdP_NE(i, j) = sum(ws_NE.*polyval(conv(dPHI_NE(i, :), dPHI_NE(j, :)), chis_NE));
        
        PP_M(i, j) = sum(ws.*polyval(conv(PHI_M(i, :), PHI_M(j, :)), chis));
        PdP_M(i, j) = sum(ws.*polyval(conv(PHI_M(i, :), dPHI_M(j, :)), chis));
        dPdP_M(i, j) = sum(ws.*polyval(conv(dPHI_M(i, :), dPHI_M(j, :)), chis));
    end
end

tol = 1e-10;
PP_NL(abs(PP_NL)<tol)=0;
PP_NE(abs(PP_NE)<tol)=0;
PP_M(abs(PP_M)<tol)=0;
PdP_NL(abs(PdP_NL)<tol)=0;
PdP_NE(abs(PdP_NE)<tol)=0;
PdP_M(abs(PdP_M)<tol)=0;
dPdP_NL(abs(dPdP_NL)<tol)=0;
dPdP_NE(abs(dPdP_NE)<tol)=0;
dPdP_M(abs(dPdP_M)<tol)=0;

figure(1)
% heatmap(PP_NL)
spy(PP_NL)
title('Lagrangian Mass Matrix')
figure(2)
% heatmap(PP_NE)
spy(PP_NE)
title('Espectral Mass Matrix')
figure(3)
% heatmap(PP_M)
spy(PP_M)
title('Modal Mass Matrix')

figure(4)
spy(PdP_NL)
title('Lagrangian Phi-dPhi Matrix')
figure(5)
spy(PdP_NE)
title('Espectral Phi-dPhi Matrix')
figure(6)
spy(PdP_M)
title('Modal Phi-dPhi Matrix')

figure(7)
spy(dPdP_NL)
title('Lagrangian Laplacian Matrix')
figure(8)
spy(dPdP_NE)
title('Espectral Laplacian Matrix')
figure(9)
spy(dPdP_M)
title('Modal Laplacian Matrix')