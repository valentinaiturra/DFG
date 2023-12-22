function fcoriolis(lat)
rad = (lat*pi)/180; %Angulo de latitud a radián
omega = 7.292*(10^(-5)); % [s^-1]
a = 6378; % [Km] radio terrestre
a = a*1000; % [m]

if (lat<=90) && (lat>= -90) 
    f = 2*omega*sin(rad);
    beta = (2*omega*cos(rad))/a;
    m = ['f = ', num2str(f),'[1/s]'];
    n = ['Beta = ' num2str(beta), '[1/ms]'];

    disp('El parámetro de coriolis está dado por ')
    disp(m)
    disp('')
    disp('El parámetro beta está dado por')
    disp(n)
else
    disp('Error, ingrese una latitud válida')
end
