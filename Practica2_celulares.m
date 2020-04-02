clear all;
clc;
%--------------------------------------------------------------------------
c11 = [1];
aux = c11;
%Potencia de generacion de chip, 2^pott
%2^5 = 32 , maxima longitud de chip 
pott = 5;
for i = 1:pott
    [m,n] = size(aux);
%     m
    ix = 0;
    %Para las filas de cada matriz
    for xi = 1:m
        cx = [aux(xi,:) ,  aux(xi,:); 
              aux(xi,:) , -aux(xi,:)];
        if(ix == 0)
            cy = cx;
            ix = 1;
        elseif(ix == 1)
            cy = [cy;cx];
            ix = 1;
        end
    end
%     cy;
    %Guardar matrices en su respectivo
    aux = cy;
    switch i
        case 1
             rmin16 = aux;
        case 2
             rmin8 = aux;
        case 3
             rmin4 = aux;
        case 4
             rmin2 = aux;
        case 5
             rmin = aux;
        otherwise
             rmin = 0;
    end
end
%--------------------------------------------------------------------------
% 
% rmin16 %c(2,1), c(2,2)
% rmin8  %c(4,1), c(4,2), c(4,3), c(4,4)
% rmin4  %c(8,1), c(8,2), c(8,3), c(8,4), c(8,5), c(8,6), c(8,7), c(8,8)...
% rmin2
% rmin

%Para generar cada usuario
%4 usuarios de tasa Rmin
rmin(29,:);
rmin(30,:);
rmin(31,:);
rmin(32,:);
%2 usuarios de tasa 2Rmin
rmin2(11,:);
rmin2(12,:);
%1 usuario de tasa 4Rmin
rmin4(4,:);

%X usuarios de tasa RminX
usuarios = [1 -1 -1  1]; %1 bits po cada uno (4 usuarios)
usuarios2= [1  1 -1 -1]; %2 bits por cada uno (2 usuarios)
usuarios4= [1 -1  1  1]; %4 bits por cada uno (1 usuario)

              %Codigo * usuario
    xc1 = rmin(29,:).*usuarios(1,1);
    xc2 = rmin(30,:).*usuarios(1,2);
    xc3 = rmin(31,:).*usuarios(1,3);
    xc4 = rmin(32,:).*usuarios(1,4);
    
    xc5x1 = rmin2(11,:).*usuarios2(1,1);
    xc5x2 = rmin2(11,:).*usuarios2(1,2);
    xc5x = [xc5x1 xc5x2];
    xc5y1 = rmin2(12,:).*usuarios2(1,3);
    xc5y2 = rmin2(12,:).*usuarios2(1,4);
    xc5y = [xc5y1 xc5y2];
    
    xc61 = rmin4(4,:).*usuarios4(1,1);
    xc62 = rmin4(4,:).*usuarios4(1,2);
    xc63 = rmin4(4,:).*usuarios4(1,3);
    xc64 = rmin4(4,:).*usuarios4(1,4);
    xc6 = [xc61 xc62 xc63 xc64];
%--------------------------------------------------------------------------    

s1d = xc1 + xc2 + xc3 + xc4 + xc5x + xc5y + xc6;
[m tn] = size(s1d);
t = 0:tn-1;

%--------------------------------------------------------------------------
%Receptor, dedispersion. Multiplicar por codigo e integrar

    des(1,:) = s1d(1,:) .* rmin(29,:);
    des(2,:) = s1d(1,:) .* rmin(30,:);
    des(3,:) = s1d(1,:) .* rmin(31,:);
    des(4,:) = s1d(1,:) .* rmin(32,:);
    
s2b = zeros(1,32);
trans = zeros(1,32);
for i = 1:4 %filas
    for ii = 1:32
            s2b(i) = s2b(i) + des(i,ii);
    end
    if  s2b(i) < 0 
        trans(i) = -1;
    elseif s2b(i) > 0 
        trans(i) = 1;
    end
end

    des2(1,:) = s1d(1:16) .* rmin2(11,:);
    des2(2,:) = s1d(17:32).* rmin2(11,:);
    des2(3,:) = s1d(1:16) .* rmin2(12,:);
    des2(4,:) = s1d(17:32).* rmin2(12,:);
    
s2b2x = zeros(1,16);
trans2x = zeros(1,16);
for i = 1:4 %filas
    for ii = 1:16
            s2b2x(i) = s2b2x(i) + des2(i,ii);
    end
    if  s2b2x(i) < 0 
        trans2x(i) = -1;
    elseif s2b2x(i) > 0 
        trans2x(i) = 1;
    end
end
    
    des4(1,:) = s1d(1:8)  .* rmin4(4,:);
    des4(2,:) = s1d(9:16) .* rmin4(4,:);
    des4(3,:) = s1d(17:24).* rmin4(4,:);
    des4(4,:) = s1d(25:32).* rmin4(4,:);

s2b4 = zeros(1,8);
trans4 = zeros(1,8);
for i = 1:4 %filas
    for ii = 1:8
            s2b4(i) = s2b4(i) + des4(i,ii);
    end
    if  s2b4(i) < 0 
        trans4(i) = -1;
    elseif s2b4(i) > 0 
        trans4(i) = 1;
    end
end
%--------------------------------------------------------------------------
%Comparativa transmisor receptor
usuarios
trans

usuarios2
trans2x

usuarios4
trans4

%--------------------------------------------------------------------------

figure(1)
stem(t,s1d);
grid on;
xlabel('t');
ylabel('amplitud');
title('s(t), señal resultante');
