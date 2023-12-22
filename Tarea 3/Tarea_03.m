clc
clear all

CTD10 = readmatrix('KN38D_P06_010.txt');
CTD12 = readmatrix('KN38D_P06_012.txt');

P10 = CTD10(:,1);
Tabla10 = CTD10(:,2);
S10 = CTD10(:,3);
lat10 = -32.490167;
lon10 = -72.512667;

P12 = CTD12(:,1);
T12 = CTD12(:,2);
S12 = CTD12(:,3);
lat12 = -32.503500;
lon12 = -72.997500;

%% Calcule la salinidad absoluta a partir de la salinidad practica dada en 
% los archivos indicados anteriormente. Use la función:  gsw_SA_from_SP

[SA10, ~] = gsw_SA_from_SP(S10,P10,lon10,lat10);
[SA12, ~] = gsw_SA_from_SP(S12,P12,lon12,lat12);

%% Calcule la temperatura conservativa (Q) a partir de q 
% (temperatura potencial). Use la función:   gsw_CT_from_pt.  Explique qué
% representa la temperatura conservativa (Q). 
% Use el documento del siguiente link:  
% http://www.teos-10.org/pubs/gsw/pdf/CT_from_pt.pdf, más otros que 
% considere útiles (incluya las referencias correspondientes). 

TP10 = gsw_pt_from_t(SA10,Tabla10,P10);
TP12 = gsw_pt_from_t(SA12,T12,P12);

CT10 = gsw_CT_from_pt(SA10,TP10);
CT12 = gsw_CT_from_pt(SA12,TP12);


% Grafique los perfiles de temperatura conservativa (versus presión) y 
% salinidad absoluta (versus presión) de las dos estaciones. Grafique la 
% temperatura de las dos estaciones en un mismo gráfico para compararlas, 
% e igualmente para la salinidad.  

% % figure(1)
% % plot(CT10,P10,'LineWidth',2)
% % hold on
% % plot(CT12,P12,'LineWidth',2)
% % grid minor
% % axis tight
% % legend('CTD10','CTD12')
% % xlabel('Temperatura Conservativa')
% % ylabel('Presión')
% % title('Perfiles de temperatura conservativa versus presión')
% % axis ij
% % 
% % figure(2)
% % plot(SA10,P10,'LineWidth',2)
% % hold on
% % plot(SA12,P12,'LineWidth',2)
% % grid minor
% % axis tight
% % legend('CTD10','CTD12')
% % xlabel('Salinidad Absoluta')
% % ylabel('Presión')
% % title('Perfiles de salinidad absoluta versus presión')
% % axis ij

%% Use la función gsw_SA_CT_plot para graficar la salinidad absoluta versus
% la temperatura conservativa. Incluya en el gráfico las isopicnas 
% correspondientes para una adecuada visualización.

%Ejemplo: gsw_SA_CT_plot(SA010, CT010,0,(25.5:0.5:27))
% Donde SA010 y CT010 son la salinidad absoluta y la temperatura 
% conservativa de la estación 010, respectivamente, 0 es la presión de 
% referencia y el último vector contiene las isopicnas a incluir en el 
% gráfico (elija usted un rango de isopicnas adecuadas para que el diagrama
% se visualice adecuadamente).


% % figure(1)
% % subplot(1,2,1)
% % gsw_SA_CT_plot(SA10, CT10,0,(19:0.5:30),'\it{S}\rm_A - {\Theta} estación CTD10')
% % subplot(1,2,2)
% % gsw_SA_CT_plot(SA12, CT12,0,(19:0.5:30),'\it{S}\rm_A - {\Theta} estación CTD12')

%% Usando la nueva ecuación de estado, calcule la densidad potencial en 
% cada estación usando como presión de referencia la superficie 
% (presión cero). 

rho_p10 = gsw_pot_rho_t_exact(SA10,Tabla10,P10,0);

rho_p12 = gsw_pot_rho_t_exact(SA12,T12,P12,0);


%Grafique los perfiles (figuras densidad potencial versus presión para 
% ambas estaciones)

% % figure()
% % plot(rho_p10,P10,'LineWidth',2)
% % hold on
% % plot(rho_p12,P12,'LineWidth',2)
% % grid minor
% % axis tight
% % legend('CTD10','CTD12')
% % xlabel('Densidad potencial')
% % ylabel('Presión')
% % title('Perfiles de densidad potencial versus presión')
% % axis ij

%% Integre numéricamente (usando trapz.m) la densidad de cada estación 
% entre la superficie y los 1000 db (note que el intervalo de la presión 
% en los datos es de 2 db)
x = [1:2:1001];
x = x';

%Buscamos las posiciones en que se encuentran la presiones hasta 1501
for i = 1:length(x)
    position(i) = find(P10 == x(i)); 
end

rho_p10 = rho_p10(position);
rho_p12 = rho_p12(position);

a = 1;
for i = 2:length(position)
    a=a+1;
    PA(a) = (P10(i-1)+P10(i))/2;
    rho_p10(a) = (rho_p10(i-1)+rho_p10(i))/2;
    rho_p12(a) = (rho_p12(i-1)+rho_p12(i))/2;
end
PA = PA'; %PA son las presiones de 0 a 1000
rho_p10 = rho_p10';
rho_p12 = rho_p12';

Q1 = trapz(PA,rho_p10);
Q2 = trapz(PA,rho_p12);

%% Obtenga la anomalía del volumen específico () en cada estación

%Usando la función: gsw_specvol_anom_standard(SA,CT, p) donde SA, 
% CT y p son salinidad absoluta, temperatura conservativa y presión
% (en decibares "db”) respectivamente.

SV10 = gsw_specvol_anom_standard(SA10,CT10,P10);

SV12 = gsw_specvol_anom_standard(SA12,CT12,P12);

%% Construya una tabla similar a la tabla de más abajo),
% usando los valores de las estaciones oceanográficas dadas

rho10 = gsw_rho(SA10,CT10,P10)-1000;
rho12 = gsw_rho(SA12,CT12,P12)-1000;

a = 1;
for i = 2:(length(P12))
    a=a+1;
    %CTD10
    datos10(a,2) = (P10(i-1)+P10(i))/2; %Presion
    datos10(a,5) = (rho10(i-1)+rho10(i))/2;
    datos10(a,6) = (SV10(i-1)+SV10(i))/2; %Volumen especifico
    datos10(a,3) = (T10(i-1)+T10(i))/2; %Temperatura
    datos10(a,4) = (S10(i-1)+S10(i))/2; %Salinidad
    %CTD12
    datos12(a,2) = (P12(i-1)+P12(i))/2; %Presion
    datos12(a,5) = (rho12(i-1)+rho12(i))/2;
    datos12(a,6) = (SV12(i-1)+SV12(i))/2; %Volumen especifico
    datos12(a,3) = (T12(i-1)+T12(i))/2; %Temperatura
    datos12(a,4) = (S12(i-1)+S12(i))/2; %Salinidad
end
%Para una presion cero, repetimos los valores de la presion uno
%CTD10
datos10(1,5) = rho10(1);
datos10(1,3) = T(1); %Temperatura
datos10(1,4) = S10(1);
datos10(1,6) = SV10(1); %Volumen especifico

%CTD12
datos12(1,5) = rho12(1);
datos12(1,3) = T12(1); %Temperatura
datos12(1,4) = S12(1);
datos12(1,6) = SV12(1); %Volumen especifico

%Buscamos donde se encuentran las presiones o profundidades estandares
depth = [0:10:30,50,76,100,126,150:50:300,400:100:1500,1750,2000:500:4500];
for i = 1:length(depth)
    d(i) = find(datos10(:,2) == depth(i));%Posiciones a esa profundidad
end

%Convierto mi matriz a una con los datos a produndidades estandar
datos10 = datos10(d,:);
datos12 = datos12(d,:);

%Determinamos el valor de la profundidad de cada dato
datos10(:,1) = -1*gsw_z_from_p(datos10(:,2),lat10);
datos12(:,1) = -1*gsw_z_from_p(datos12(:,2),lat12);


a = 1;
for i = 2:length(datos10(:,1))
    a = a+1;
    %Determinamos el volumen especifico promedio
    datos10(a,7) = (datos10((i-1),6) + datos10(i,6))/2;
    datos12(a,7) = (datos12((i-1),6) + datos12(i,6))/2;
    %Determinamos diferencia de presion, es igual para CTD10 y CTD12
    DP(a) = (datos10(i,2) - datos10((i-1),2))*10000; %Difernecia de presion en Pascales
end

datos10(1,7) = NaN;
datos12(1,7) = NaN;
DP(1) = NaN;
format shortG
for i = 1:length(datos10(:,1))
    %La ultima columna de la tabla es volumen especifico promedio y
    %diferencia de presion
    datos10(i,8) = datos10(i,7)*DP(i);
    datos12(i,8) = datos12(i,7)*DP(i);
end
Tabla10 = table(datos10(:,1),datos10(:,2),datos10(:,3),datos10(:,4),datos10(:,5),datos10(:,6),datos10(:,7),datos10(:,8),'VariableNames',["depth [m]","presion [dbar]","Temperatura [°C]","Salinidad","Densidad [kg m^-3]","Delta [m^3/kg]","Delta prom [m^3/kg]","Multiplicacion [m^3*Pa/kg]"]);

Tabla12 = table(datos12(:,1),datos12(:,2),datos12(:,3),datos12(:,4),datos12(:,5),datos12(:,6),datos12(:,7),datos12(:,8),'VariableNames',["depth [m]","presion [dbar]","Temperatura [°C]","Salinidad","Densidad [kg m^-3]","Delta [m^3/kg]","Delta prom [m^3/kg]","Multiplicacion [m^3*Pa/kg]"]);

writetable(Tabla10, 'tabla10.xlsx');
writetable(Tabla12, 'tabla12.xlsx');

%% Integre numéricamente las anomalías del volumen específico entre los distintos pares de isóbaras
% (entre 1500 db y cada nivel de presión hasta la superficie) para obtener la anomalía geopotencial (puede usar
% los valores de la tabla directamente sumando la penúltima columna, como en la tabla del ejemplo. Se debe
% tener cuidado en la integración si las presiones no están equi-espaciadas). Si es necesario extrapole
% linealmente los valores hasta p = 0.
 
AGP10 = datos10(2:24,8);
AGP12 = datos12(2:24,8);

AGP10 = flipud(AGP10);
AGP12 = flipud(AGP12);

AGP10 = cumsum(AGP10);
AGP12 = cumsum(AGP12);

AGP10 = flipud(AGP10);
AGP12 = flipud(AGP12);

AGP10(24) = 0;
AGP12(24) = 0;

%% Dadas las posiciones geográficas calcule la distancia entre las estaciones. Indique y explique el método
% que utilizó para calcular distancia entre dos puntos sobre la Tierra. Note que el paquete The Gibbs SeaWater
% (GSW) Oceanographic Toolbox tiene una opción.

%Utilizando teor de pitagoras
LONG = ((lon10 -lon12)*96486).^2; %[metros]
LAT = ((lat10 - lat12)*110852).^2; %[metros]
dis = sqrt((LONG) + (LAT));

%Confirmando con la funcion
long = [lon10 lon12];
lat = [lat10 lat10];

distance = gsw_distance(long,lat);
%% Obtenga la variación de la velocidad geostrófica entre dos niveles de 
% presión p1 y p2 mediante larelación Donde los subíndices A y B denotan 
% las dos estaciones separadas una distancia L y f es el parámetro de 
% Coriolis evaluado a una latitud promedio entre ambas estaciones. 
% Suponiendoque la velocidad en 1500 db es cero, obtenga la velocidad V 
% en los distintos niveles e incluya estos valores en una Tabla.

%Parametros necesarios
latprom = ((lat12)+(lat10))/2;
rad = (latprom*pi)/180; %Angulo de latitud a radián
omega = 7.292*(10^(-5)); % [s^-1]
f = 2*omega*sin(rad); %parametro de coriolis
L = distance;

%Determinamos primeramente la variacion de velocidad geostrofica 
for i = 1:(length(AGP12))
    dV(i,1) = (1/(L*f))*(AGP10(i) - AGP12(i)) 
end
dV = flipud(dV)

V(1) = 0; %Sabemos que en 1500 dbar la velocidad es cero
for i = 2:(length(dV))
    V(i,1) = dV(i) + V(i-1)
end
T= table(datos10(1:24,2),AGP10,flipud(V(1:end)));

%% Grafique la velocidad geostrófica entre las estaciones A y B versus 
% presión. INTERPRETE Y COMENTE EL GRÁFICO.

plot(flipud(V(1:end)),datos10(1:24,2))
axis ij


