clc
clear figures

N = 5;
dx = .01;

equis = -1:dx:1';
Hp = zeros(length(-1:dx:1), N+1);
Hps = Base_NodalEspectral(N);
for n = 0:N
    Hp(:,n+1) = polyval(Hps(n+1,:),equis);
end

figure(1)
plot(equis, Hp)%, equis, PJAC2)
legend(string(0:N),'Interpreter','LaTex','Fontsize',12)%,'location','best')
xlabel('$\xi$','Interpreter','LaTex','Fontsize',12)
ylabel('$h_p(\xi)$','Interpreter','LaTex','Fontsize',12)
%title('Aproximación Nodal Espectral','Interpreter','LaTex','Fontsize',12)

Lp = zeros(length(-1:dx:1), N+1);
Lps = Base_NodalLagrange(N);
for n = 0:N
    Lp(:,n+1) = polyval(Lps(n+1,:),equis);
end

figure(2)
plot(equis, Lp)%, equis, PJAC2)
legend(string(0:N),'Interpreter','LaTex','Fontsize',12)%,'location','best')
xlabel('$\xi$','Interpreter','LaTex','Fontsize',12)
ylabel('$\Phi_p(\xi)$','Interpreter','LaTex','Fontsize',12)
%title('Aproximación Nodal Lagrange','Interpreter','LaTex','Fontsize',12)


Pp = zeros(length(-1:dx:1), N+1);
Pps = Base_Modal(N, 1, 1);
for n = 0:N
    Pp(:,n+1) = polyval(Pps(n+1,:),equis);
end

figure(3)
plot(equis, Pp)%, equis, PJAC2)
legend(string(0:N),'Interpreter','LaTex','Fontsize',12)%,'location','best')
xlabel('$\xi$','Interpreter','LaTex','Fontsize',12)
ylabel('$\psi_p(\xi)$','Interpreter','LaTex','Fontsize',12)
%title('Aproximación Modal Espectral','Interpreter','LaTex','Fontsize',12)
