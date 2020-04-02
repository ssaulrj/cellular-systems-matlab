clc; clear all; close all;
r = 2000; %radio
%Puntos de origen----------------------------------------------------------
xs = [0    0 3000  3000     0 -3000 -3000];
ys = [0 3465 1750 -1750 -3465 -1750 1750];
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
plot(0,0,'xb',0,3465,'xb',3000,1750,'xb',3000,-1750,'xb',0,-3465,'xb',-3000,-1750,'xb',-3000,1750,'xb'); %centros
% axis square;
grid on;
hold on;

%Distancia euclideana
% sqrt( ((x-x2).^2) + ((y-y2).^2) )
% Para generar cada usuario

for i = 1:100
%     rx = floor((rand()*(2*r+1)) - 2000);
%     ry = floor((rand()*(2*r+1)) - 2000);
    rx = floor((rand()*(2*1750+1)) - 1750);
    ry = floor((rand()*(2*1750+1)) - 1750);
    ax(1,i) = rx;
    ay(1,i) = ry;
    dxy(1,i) = floor( sqrt( ((xs(1)-ax(1,i))^2) + ((ys(1)-ay(1,i))^2) ) );
    plot(ax(1,i),ay(1,i),'o');
    hold on;
end
ax
ay
dxy
xlim([-6000 6000])
ylim([-6000 6000])
%Determinar la probabilidad de Outage--------------------------------------
% Poutage = (numsinservicio)/(totalusuarios);
