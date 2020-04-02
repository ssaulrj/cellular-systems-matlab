clc; clear all; close all;
r = 6000; %radio
%Puntos de origen----------------------------------------------------------
xs = [0    0 r*(3/2)  r*(3/2)     0 -r*(3/2) -r*(3/2)];
ys = [0 r*1.7325 r*(7/8) -r*(7/8) -r*1.7325 -r*(7/8) r*(7/8)];
                     %N,r,cx,cy
[x,y]  = hexacelulares(6,r,xs(1),ys(1));
[x2,y2]= hexacelulares(6,r,xs(2),ys(2));
[x3,y3]= hexacelulares(6,r,xs(3),ys(3));
[x4,y4]= hexacelulares(6,r,xs(4),ys(4));
[x5,y5]= hexacelulares(6,r,xs(5),ys(5));
[x6,y6]= hexacelulares(6,r,xs(6),ys(6));
[x7,y7]= hexacelulares(6,r,xs(7),ys(7));

figure(1)
plot(x,y,'k',x2,y2,'k',x3,y3,'k',x4,y4,'k',x5,y5,'k',x6,y6,'k',x7,y7,'k');
grid on;
hold on;
plot(xs(1),ys(1),'xb',xs(2),ys(2),'xb',xs(3),ys(3),'xb',xs(4),ys(4),'xb',xs(5),ys(5),'xb',xs(6),ys(6),'xb',xs(7),ys(7),'xb'); %centros
% axis square;
grid on;
hold on;

%Distancia euclideana
% sqrt( ((x-x2).^2) + ((y-y2).^2) )
% Para generar cada usuario
Ptx = 2:2:10;

totalusuarios = 1000;
% ensombrecimiento = 5:8;
enx = 5:8;
for j = 5:8 %ensombrecimiento
    numsinservicio = 0;
    for i = 1:totalusuarios
    %     rx = floor((rand()*(2*r+1)) - 2000);
    %     ry = floor((rand()*(2*r+1)) - 2000);
        while 1
            rx = floor((rand()*(2*r*(7/8)+1)) - r*(7/8));
            ry = floor((rand()*(2*r*(7/8)+1)) - r*(7/8));
        %     if(ry>)
            dist = floor( sqrt( ((xs(1)-rx)^2) + ((ys(1)-ry)^2) ) );
            if dist > r*(7/8)

            else 
                break;
            end
        end
        ax(1,i) = rx;
        ay(1,i) = ry;
    %     dxy0(1,i) = floor( sqrt( ((xs(1)-ax(1,i))^2) + ((ys(1)-ay(1,i))^2) ) ); % distancia de punto 0 a coordenadas
        dxy0(1,i) = sqrt( ((xs(1)-ax(1,i))^2) + ((ys(1)-ay(1,i))^2) ); % distancia de punto 0 a coordenadas
        pxt0(1,i) = 40 + 12 + 2 - (10*4*log10( dxy0(1,i) )) - j; %potencia rx
        if pxt0(1,i) < -108
            numsinservicio = numsinservicio + 1;
        end 
    %     dxy1(1,i) = floor( sqrt( ((xs(2)-ax(1,i))^2) + ((ys(2)-ay(1,i))^2) ) );
    %     dxy2(1,i) = floor( sqrt( ((xs(3)-ax(1,i))^2) + ((ys(3)-ay(1,i))^2) ) );
    %     dxy3(1,i) = floor( sqrt( ((xs(4)-ax(1,i))^2) + ((ys(4)-ay(1,i))^2) ) );
    %     dxy4(1,i) = floor( sqrt( ((xs(5)-ax(1,i))^2) + ((ys(5)-ay(1,i))^2) ) );
    %     dxy5(1,i) = floor( sqrt( ((xs(6)-ax(1,i))^2) + ((ys(6)-ay(1,i))^2) ) );
    %     dxy6(1,i) = floor( sqrt( ((xs(7)-ax(1,i))^2) + ((ys(7)-ay(1,i))^2) ) );
        plot(ax(1,i),ay(1,i),'o');
        hold on;
    end
    Poutage = (numsinservicio)/(totalusuarios);
    ensombrecimiento(1,j-4) = Poutage;
    xlim([-r*3 r*3])
    ylim([-r*3 r*3])
end 

numsinservicio;

% figure(2)
% plot(enx,ensombrecimiento);
