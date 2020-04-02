clc; clear all; close all;
r = 5000; %radio

%Puntos de origen----------------------------------------------------------
xs = [0    0 r*(3/2)  r*(3/2)     0 -r*(3/2) -r*(3/2)];
ys = [0 r*1.7325 r*(7/8) -r*(7/8) -r*1.7325 -r*(7/8) r*(7/8)];
                      %N,r, cx  , cy
[x,y]  = hexacelulares(6,r,xs(1),ys(1));
[x2,y2]  = hexacelulares(6,r,xs(2),ys(2));
[x3,y3]  = hexacelulares(6,r,xs(3),ys(3));
[x4,y4]  = hexacelulares(6,r,xs(4),ys(4));
[x5,y5]  = hexacelulares(6,r,xs(5),ys(5));
[x6,y6]  = hexacelulares(6,r,xs(6),ys(6));
[x7,y7]  = hexacelulares(6,r,xs(7),ys(7));
figure(1)
plot(x,y,'k',x2,y2,'k',x3,y3,'k',x4,y4,'k',x5,y5,'k',x6,y6,'k',x7,y7,'k');
grid on;
hold on;
plot(xs(1),ys(1),'xb',xs(2),ys(2),'xb',xs(3),ys(3),'xb',xs(4),ys(4),'xb',xs(5),ys(5),'xb',xs(6),ys(6),'xb',xs(7),ys(7),'xb'); %centros
grid on;
hold on;
xq = (4500*randn(1,2500));
yq = (4500*randn(1,2500));
[in,on] = inpolygon(xq,yq,x,y);
xq2 = (4500*randn(1,2500));
yq2 = (4500*randn(1,2500));
[in2,on2] = inpolygon(xq2,yq2,x2,y2);
xq3 = (4500*randn(1,2500));
yq3 = (4500*randn(1,2500));
[in3,on3] = inpolygon(xq3,yq3,x3,y3);
xq4 = (4500*randn(1,2500));
yq4 = (4500*randn(1,2500));
[in4,on4] = inpolygon(xq4,yq4,x4,y4);
xq5 = (4500*randn(1,2500));
yq5 = (4500*randn(1,2500));
[in5,on5] = inpolygon(xq5,yq5,x5,y5);
xq6 = (4500*randn(1,2500));
yq6 = (4500*randn(1,2500));
[in6,on6] = inpolygon(xq6,yq6,x6,y6);
xq7 = (4500*randn(1,2500));
yq7 = (4500*randn(1,2500));
[in7,on7] = inpolygon(xq7,yq7,x7,y7);
numinside1=xq(in);
numinside2=xq(in2);
numinside3=xq(in3);
numinside4=xq(in4);
numinside5=xq(in5);
numinside6=xq(in6);
numinside7=xq(in7);
hold on;
plot(xq(in),yq(in),'mo',xq(in2),yq(in2),'co',xq(in3),yq(in3),'ro',xq(in4),yq(in4),'go',xq(in5),yq(in5),'bo',xq(in6),yq(in6),'ko',xq(in7),yq(in7),'yo') % points inside

todosx = [xq(in) xq(in2) xq(in3) xq(in4) xq(in5) xq(in6) xq(in7)];
todosy = [yq(in) yq(in2) yq(in3) yq(in4) yq(in5) yq(in6) yq(in7)];
figure(2)
plot(todosx, todosy,'x');
hold on;
plot(x,y,'k',x2,y2,'k',x3,y3,'k',x4,y4,'k',x5,y5,'k',x6,y6,'k',x7,y7,'k');
hold on; 

Pottx= 40;
Gtx = 12;
Grx = 2;
ensombrecimiento = 7;
alfa = 4;

usuariosTotales = length(todosx);
for centrosiu = 1:7
    for iu = 1:usuariosTotales
        dxy(centrosiu,iu) = floor( sqrt( ((xs(centrosiu)-todosx(iu))^2) + ((ys(centrosiu)-todosy(iu))^2) ) );
        pxt(centrosiu,iu) = Pottx + Gtx + Grx - (10*4*log10( dxy(centrosiu,iu) )) - ensombrecimiento; %potencia rx
    end
end

for iu = 1:usuariosTotales
        pxtmax = max(pxt);
end

[M , I] = max(pxt);

for iu = 1:usuariosTotales
        if I(iu) == 1
            plot(todosx(iu),todosy(iu),'mo');
            hold on;
            usuarioscentralesx(iu) = todosx(iu);
            usuarioscentralesy(iu) = todosy(iu);
        end
end
xlim([-r*8 r*8])
ylim([-r*8 r*8])
hold on;

for ks = 1:4
    figure(3)
    %factores K de reuso 1, 3, 7 y 13
    iii = [1 1 2 3];
    j = [0 1 1 1]
    %1  i = 1, j = 0
    %3 i = 1, j = 1
    %7 i = 2, j = 1
    %13 i = 3, j = 1
    %     iii = 3;
    %     j = 0;
    K = (iii(ks)^2)+(iii(ks)*j(ks))+(j(ks)^2);
    apotema = sqrt( (r^2) - ((r/2)^2) );

    xsx = 0;
    ysy = 0;
    xsx2 = 0;
    ysy2 = 0;
    xsx3 = 0;
    ysy3 = 0;
    xsx4 = 0;
    ysy4 = 0;
    xsx5 = 0;
    ysy5 = 0;
    xsx6 = 0;
    ysy6 = 0;

    aladerecha = r*1.5;
    aarriba = apotema;
    for ix = 1:iii(ks)+1
        plot(xsx,ysy,'k*');
        hold on;
        enx = xsx;
        eny = ysy;
        xsx = xsx+aladerecha;
        ysy = ysy+aarriba;

        plot(xsx2,ysy2,'k*');
        hold on;
        enx2 = xsx2;
        eny2 = ysy2;
        xsx2 = xsx2+aladerecha;
        ysy2 = ysy2-aarriba;

        plot(xsx3,ysy3,'k*');
        hold on;
        enx3 = xsx3;
        eny3 = ysy3;
        xsx3 = xsx3-aladerecha;
        ysy3 = ysy3+aarriba;

        plot(xsx4,ysy4,'k*');
        hold on;
        enx4 = xsx4;
        eny4 = ysy4;
        xsx4 = xsx4-aladerecha;
        ysy4 = ysy4-aarriba;

        plot(xsx5,ysy5,'k*');
        hold on;
        enx5 = xsx5;
        eny5 = ysy5;
        ysy5 = ysy5+(apotema*2);

        plot(xsx6,ysy6,'k*');
        hold on;
        enx6 = xsx6;
        eny6 = ysy6;
        ysy6 = ysy6-(apotema*2);
    end
    hold on;

    %anillo de interferencia
    interx = [enx enx2+(j(ks)*apotema) enx6+(j(ks)*r*1.5) enx4 enx3-(j(ks)*r*1.5) enx5-(j(ks)*r*1.5)];
    intery = [eny+(j(ks)*(apotema*2)) eny2+(j(ks)*r*1.5) eny6-(j(ks)*apotema) eny4-(j(ks)*(apotema*2))  eny3-(j(ks)*apotema)  eny5+(j(ks)*apotema)];
    plot(interx,intery,'r*')
    hold on;
    plot(interx, intery, enx, eny+(j(ks)*(apotema*2)))
    hold on;

    tamcenu = length(usuarioscentralesx);
    for iu = 1:tamcenu
            dxypot(1,iu) = floor( sqrt( ((xs(1)-usuarioscentralesx(iu))^2) + ((ys(1)-usuarioscentralesy(iu))^2) ) );
            potsenial(1,iu) = (Pottx * Gtx * Grx)/((dxypot(1,iu)^alfa)*(10^(ensombrecimiento/10))); %potencia 1
    end
    potsenialx(1,:) = zeros(1,tamcenu);
    for centrosiu = 1:6
        for iu = 1:tamcenu
            dxypot(centrosiu+1,iu) = floor( sqrt( ((interx(centrosiu)-usuarioscentralesx(iu))^2) + ((intery(centrosiu)-usuarioscentralesy(iu))^2) ) );

            potsenial(centrosiu+1,iu) = (Pottx * Gtx * Grx)/((dxypot(centrosiu+1,iu)^alfa)*(10^(ensombrecimiento/10))); %potencia demas
        end
        potsenialx(1,:) = potsenialx(1,:) + potsenial (centrosiu+1,:);
    end

    SIR = potsenial(1,:) ./ potsenialx(1,:);
    SIRdb = 10*log10(SIR);
    SIRdb20 = SIRdb(1:20);

    figure(4)
    subplot(2,2,ks);
    hist(SIRdb20);
    title(ks);
    grid on;
end
% 
% Puntosigx
% for ii = 1:length(c_x)
%         dxy0(1,ii) = sqrt( ((xs(1)-c_x(1,ii))^2) + ((ys(1)-c_y(1,ii))^2) ); % distancia de punto 0 a coordenadas
%         pxt0(1,ii) = Pottx + Gtx + Grx - (10*4*log10( dxy0(1,ii) )) - j; %potencia rx
%         if pxt0(1,ii) < -108
%            numsinservicio = numsinservicio + 1;
%         end 
%         aux(a,:)=numsinservicio;
% end
    
% % Para generar cada usuario
% % Ptx = 2:2:10;
% Ptx = -[33 36 39 40];
% Pottx=40;
% Gtx = 7;
% Grx=2;
% ro = 7:2:13;
% % ensombrecimiento = 5:8;
% enx = 5:8;
% 
% j=7;
% % for j = 5:8 %ensombrecimiento
%     numsinservicio = 0;
%     totalusuarios = 100;
%     for i = 1:7
%         c_x = r-rand(1, 3*totalusuarios)*2*r;
%         c_y = r-rand(1, 3*totalusuarios)*2*r;
% %         plot(c_x,c_y,'x');
%         [x,y]  = hexacelulares(6,r,xs(i),ys(i));
%         figure(1)
%         plot(x,y,'k'); %contorno
%         hold on;
%         plot(xs(i),ys(i),'xb'); %centros
%         hold on;
%         IN = inpolygon(c_x, c_y, x,y);
% %         plot(c_x,c_y,'x');
%         %1 inside, 0
%         %drop nodes outside the hexagon
%         c_x = c_x(IN);
%         c_y = c_y(IN);
%         %choose only N points
%         idx = randperm(length(c_x))
%         c_x = c_x(idx(1:totalusuarios));
%         c_y = c_y(idx(1:totalusuarios));
%         plot(c_x,c_y,'o');
%         hold on;
%         %    dxy0(1,i) = floor( sqrt( ((xs(1)-ax(1,i))^2) + ((ys(1)-ay(1,i))^2) ) ); % distancia de punto 0 a coordenadas
%         for ii = 1:length(c_x)
%             dxy0(1,ii) = sqrt( ((xs(i)-c_x(1,ii))^2) + ((ys(i)-c_y(1,ii))^2) ); % distancia de punto 0 a coordenadas
%             pxt0(1,ii) = Pottx + Gtx + Grx - (10*4*log10( dxy0(1,ii) )) - j; %potencia rx
%             if pxt0(1,ii) < -108
%                 numsinservicio = numsinservicio + 1;
%             end 
%         end
%     %     dxy1(1,i) = floor( sqrt( ((xs(2)-ax(1,i))^2) + ((ys(2)-ay(1,i))^2) ) );
%     %     dxy2(1,i) = floor( sqrt( ((xs(3)-ax(1,i))^2) + ((ys(3)-ay(1,i))^2) ) );
%     %     dxy3(1,i) = floor( sqrt( ((xs(4)-ax(1,i))^2) + ((ys(4)-ay(1,i))^2) ) );
%     %     dxy4(1,i) = floor( sqrt( ((xs(5)-ax(1,i))^2) + ((ys(5)-ay(1,i))^2) ) );
%     %     dxy5(1,i) = floor( sqrt( ((xs(6)-ax(1,i))^2) + ((ys(6)-ay(1,i))^2) ) );
%     %     dxy6(1,i) = floor( sqrt( ((xs(7)-ax(1,i))^2) + ((ys(7)-ay(1,i))^2) ) );
%     end
%     Poutage = (numsinservicio)/(totalusuarios);
% %     ensombrecimiento(1,j-4) = Poutage;
%     xlim([-r*3 r*3])
%     ylim([-r*3 r*3])
% % end 
% 
% numsinservicio;
% % 
% % figure(2)
% % plot(enx,ensombrecimiento);
% % xlabel('\wp');
% % ylabel('P_out');
% % grid on;
% % 
% % figure(3);
% % plot(Ptx,ensombrecimiento);
% % grid on;
% % xlabel('P_tx');
% % ylabel('P_out');
