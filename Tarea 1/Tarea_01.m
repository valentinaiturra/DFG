% Ejercicio 01 Dinámica de Fluidos Geofísicos
%
% Parámetros
g = 9.82; % [m/s^2] gravedad
alfa = 15; %[grados] inclinación del plano en grados
alfa = alfa*pi/180; %[rad] alfa en radianes

% Frecuencia angular (de rotación) del plano inclinado
omega = 0.5; % [s-1]
f = 2*omega*cos(alfa); %[s^-1] 2 * componente vertical de omega
ge = g*sin(alfa); %[m/s^2], 2 * componente de g en la dir del plano
%% Calcula y grafica las velocidades y trayectorias

t= (0:0.1:10)';
u = (ge/f)*cos(f*t)-ge/f; %[m/s] expresión para la componente x de la velocidad
v = -(ge/f)*sin(f*t); %expresión para la componente x de la velocidad

% Trayectorias: posición x, y en función del tiempo
x = (ge*sin(f*t))/(f^2)-(ge*t/f); %expresión para la componente x de la velocidad
y = (ge*cos(f*t))/(f^2)-(ge/f); %expresión para la componente x de la velocidad

% Trayectorias calculadas integrando numéricamente las velocidades
xn = cumtrapz(t,u); % cumtrapz es la función para integrar numeric.
yn = cumtrapz(t,v);

% Grafica velocidades
figure
plot(t,u,'.-b')
hold on
plot(t,v,'.-r')
axis tight
legend('componente horizontal','Componente vertical',Location='northeast')
xlabel('Tiempo [s]')
ylabel('Velocidad [m/s]')
title('Velocidad de la partícula')

% Grafica trayectoria
figure
plot(x,y,'b.-')
hold on
plot(xn,yn,'r.-')
axis tight
legend('Integrando manualmente','Integrando con regla trapezoidal',Location='northeast')
xlabel('Trayectoria en eje x [m]')
ylabel('Trayectoria en eje y [m]')
title('Trayectoria de la partícula')