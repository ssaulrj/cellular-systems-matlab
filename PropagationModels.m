clear all;
clc;
close all;
%--------------------------------------------------------------------------
d = [20:0.05:100]; %m
f = 800000000; %mhz
Gt = 1;
Gr = 1;
Ht = 10;
Hr = 10;
Pt = 1; %W
c = 3.*10^8;
landa = c/f;

% Perdidas de propagacion 
% Lespaciolibre = -10*log(Gt) - 10*log(Gr) + 20*log(f) + 20*log(d) - 147.56;
% Lsuperficiereflejada = -10*log(Gt) - 10*log(Gr) - 20*log(Ht) - 20*log(Hr) + 40*log(d);
Pespaciolibre = 10*log(Pt) + 10*log(Gt) + 10*log(Gr) + 20*log((c)./(4*pi*f*d));
% Psuperficiereflejada = 10*log(Gt) + 10*log(Gr) + 20*log((Ht.*Hr)./(d.^2));
Psuperficiereflejada = 10*log(4) + 10*log(Pt) + 20*log((landa)./(4*pi*d)) + 10*log(Gt) + 10*log(Gr) + 20*log(sin((2*pi*Ht*Hr)./(landa*d)));

figure(1);
subplot(2,1,1);
plot(d, Pespaciolibre);
xlabel('d [m]');
ylabel('[dB]');
title('(1) Potencia recibida, Espacio libre, ecuación 1');
subplot(2,1,2);
plot(d, Psuperficiereflejada);
xlabel('d [m]');
ylabel('[dB]');
title('(1) Potencia recibida, Superficie reflejante, ecuación 2');

%--------------------------------------------------------------------------
d2 = [100:100:3000]; %m
Pespaciolibre2 = 10*log(Pt) + 10*log(Gt) + 10*log(Gr) + 20*log((c)./(4*pi*f*d2));
Pecuacion2 = 10*log(Pt) + 10*log(Gt) + 10*log(Gr)  + 20*log((Ht.*Hr)./(d2.^2)); %Si d >> Ht, Hr

figure(2);
subplot(2,1,1);
plot(d2, Pespaciolibre2);
xlabel('d [m]');
ylabel('[dB]');
title('(2) Potencia recibida, Espacio libre, ecuación 1');
subplot(2,1,2);
plot(d2, Pecuacion2);
xlabel('d [m]');
ylabel('[dB]');
title('(2) Potencia recibida, Superficie reflejante, ecuación 3');

%--------------------------------------------------------------------------
Ht3 = 50;
Hr3 = 50;
Pespaciolibre3 = 10*log(Pt) + 10*log(Gt) + 10*log(Gr) + 20*log((c)./(4*pi*f*d2));
Pecuacion3 = 10*log(Pt) + 10*log(Gt) + 10*log(Gr)  + 20*log((Ht3.*Hr3)./(d2.^2));

figure(3);
subplot(2,1,1);
plot(d2, Pespaciolibre3);
xlabel('d [m]');
ylabel('[dB]');
title('(3) Potencia recibida, Espacio libre, ecuación 1, H_T = H_R = 50 m');
subplot(2,1,2);
plot(d2, Pecuacion3);
xlabel('d [m]');
ylabel('[dB]');
title('(3) Potencia recibida, Superficie reflejante, ecuación 3, H_T = H_R = 50 m');
