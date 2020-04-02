clc; clear all; close all;
r = 3000; %radio

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
% figure(1)
% plot(x,y,'k',x2,y2,'k',x3,y3,'k',x4,y4,'k',x5,y5,'k',x6,y6,'k',x7,y7,'k');
% grid on;
% hold on;
% plot(xs(1),ys(1),'xb',xs(2),ys(2),'xb',xs(3),ys(3),'xb',xs(4),ys(4),'xb',xs(5),ys(5),'xb',xs(6),ys(6),'xb',xs(7),ys(7),'xb'); %centros
% grid on;
% hold on;
usuarioscelda = 500;
xq = (4500*randn(1,usuarioscelda));
yq = (4500*randn(1,usuarioscelda));
[in,on] = inpolygon(xq,yq,x,y);
xq2 = (4500*randn(1,usuarioscelda));
yq2 = (4500*randn(1,usuarioscelda));
[in2,on2] = inpolygon(xq2,yq2,x2,y2);
xq3 = (4500*randn(1,usuarioscelda));
yq3 = (4500*randn(1,usuarioscelda));
[in3,on3] = inpolygon(xq3,yq3,x3,y3);
xq4 = (4500*randn(1,usuarioscelda));
yq4 = (4500*randn(1,usuarioscelda));
[in4,on4] = inpolygon(xq4,yq4,x4,y4);
xq5 = (4500*randn(1,usuarioscelda));
yq5 = (4500*randn(1,usuarioscelda));
[in5,on5] = inpolygon(xq5,yq5,x5,y5);
xq6 = (4500*randn(1,usuarioscelda));
yq6 = (4500*randn(1,usuarioscelda));
[in6,on6] = inpolygon(xq6,yq6,x6,y6);
xq7 = (4500*randn(1,usuarioscelda));
yq7 = (4500*randn(1,usuarioscelda));
[in7,on7] = inpolygon(xq7,yq7,x7,y7);
numinside1=xq(in);
numinside2=xq(in2);
numinside3=xq(in3);
numinside4=xq(in4);
numinside5=xq(in5);
numinside6=xq(in6);
numinside7=xq(in7);
% hold on;
% plot(xq(in),yq(in),'mo',xq2(in2),yq2(in2),'co',xq3(in3),yq3(in3),'ro',xq4(in4),yq4(in4),'go',xq5(in5),yq5(in5),'bo',xq6(in6),yq6(in6),'ko',xq7(in7),yq7(in7),'yo') % points inside

todosx = [xq(in) xq2(in2) xq3(in3) xq4(in4) xq5(in5) xq6(in6) xq7(in7)];
todosy = [yq(in) yq2(in2) yq3(in3) yq4(in4) yq5(in5) yq6(in6) yq7(in7)];
% figure(2)
% plot(todosx, todosy,'x');
% hold on;
% plot(x,y,'k',x2,y2,'k',x3,y3,'k',x4,y4,'k',x5,y5,'k',x6,y6,'k',x7,y7,'k');
% hold on; 

Pottx= 40;
Gtx = 15.84;
Grx = 1.58;
ensombrecimiento1 = 7;
ensombrecimiento = 7;
ensombrecimientoint = 0:5;
alfa = 4;

usuariosTotales = length(todosx);
for centrosiu = 1:7
    for iu = 1:usuariosTotales
        dxy(centrosiu,iu) = floor( sqrt( ((xs(centrosiu)-todosx(iu))^2) + ((ys(centrosiu)-todosy(iu))^2) ) );
        pxt(centrosiu,iu) = Pottx + Gtx + Grx - (10*4*log10( dxy(centrosiu,iu) )) - ensombrecimiento1; %potencia rx
    end
end

for iu = 1:usuariosTotales
        pxtmax = max(pxt);
end

[M , I] = max(pxt);


for iu = 1:usuariosTotales
        if I(iu) == 1
%             plot(todosx(iu),todosy(iu),'mo');
%             hold on;
            usuarioscentralesx(iu) = todosx(iu);
            usuarioscentralesy(iu) = todosy(iu);
        end
end
% xlim([-r*1.5 r*1.5])
% ylim([-r*1.5 r*1.5])
% hold on;

for ks = 1:1
%     figure(3)
%     plot(x,y,'k',x2,y2,'k',x3,y3,'k',x4,y4,'k',x5,y5,'k',x6,y6,'k',x7,y7,'k');
%     hold on; 
    %factores K de reuso 1, 3, 7 y 13
    iii = [1 1 2 3];
    j = [0 1 1 1];
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
%         plot(xsx,ysy,'k*');
%         hold on;
        enx = xsx;
        eny = ysy;
        xsx = xsx+aladerecha;
        ysy = ysy+aarriba;

%         plot(xsx2,ysy2,'k*');
%         hold on;
        enx2 = xsx2;
        eny2 = ysy2;
        xsx2 = xsx2+aladerecha;
        ysy2 = ysy2-aarriba;

%         plot(xsx3,ysy3,'k*');
%         hold on;
        enx3 = xsx3;
        eny3 = ysy3;
        xsx3 = xsx3-aladerecha;
        ysy3 = ysy3+aarriba;

%         plot(xsx4,ysy4,'k*');
%         hold on;
        enx4 = xsx4;
        eny4 = ysy4;
        xsx4 = xsx4-aladerecha;
        ysy4 = ysy4-aarriba;

%         plot(xsx5,ysy5,'k*');
%         hold on;
        enx5 = xsx5;
        eny5 = ysy5;
        ysy5 = ysy5+(apotema*2);

%         plot(xsx6,ysy6,'k*');
%         hold on;
        enx6 = xsx6;
        eny6 = ysy6;
        ysy6 = ysy6-(apotema*2);
    end
%     hold on;

    %anillo de interferencia
    interx = [enx enx2+(j(ks)*r*1.5) enx6+(j(ks)*r*1.5) enx4 enx3-(j(ks)*r*1.5) enx5-(j(ks)*r*1.5) enx];
    intery = [eny+(j(ks)*(apotema*2)) eny2+(j(ks)*apotema) eny6-(j(ks)*apotema) eny4-(j(ks)*(apotema*2))  eny3-(j(ks)*apotema)  eny5+(j(ks)*apotema) eny+(j(ks)*(apotema*2))];
%     plot(interx,intery,'r*');
%     hold on;
%     plot(interx, intery);
%     hold on;

    tamcenu = length(usuarioscentralesx);
    for iu = 1:tamcenu
            dxypot(1,iu) = floor( sqrt( ((xs(1)-usuarioscentralesx(iu))^2) + ((ys(1)-usuarioscentralesy(iu))^2) ) );
            potsenial(1,iu) = (Pottx * Gtx * Grx)./((dxypot(1,iu)^alfa)*(10^(ensombrecimiento/10))); %potencia 1
    end
    potsenialx(1,:) = zeros(1,tamcenu);
    for centrosiu = 1:6
        for iu = 1:tamcenu
            dxypot(centrosiu+1,iu) = floor( sqrt( ((interx(centrosiu)-usuarioscentralesx(iu))^2) + ((intery(centrosiu)-usuarioscentralesy(iu))^2) ) );

            potsenial(centrosiu+1,iu) = (Pottx * Gtx * Grx)./((dxypot(centrosiu+1,iu)^alfa)*(10^(ensombrecimientoint(centrosiu)/10))); %potencia demas
        end
        potsenialx(1,:) = potsenialx(1,:) + potsenial (centrosiu+1,:);
    end

    SIR = potsenial(1,:) ./ potsenialx(1,:);
    SIRdb = 10.*log10(SIR);
    

%     figure(4)
% %     subplot(2,2,ks);
%     hist(SIRdb,40);
%     title('Factor de reuso k =1');
%     grid on;
end

numts = zeros(1,16);

for ts = 1:tamcenu
    if SIRdb(ts) < -7.28
        SIRts(ts) = 0;
        numts(1) = numts(1) + 1;
    elseif SIRdb(ts) > -7.28 && SIRdb(ts) < -4.78
        SIRts(ts) = 19.1898;
        numts(2) = numts(2) + 1;
    elseif SIRdb(ts) > -4.78 && SIRdb(ts) < -2.04
        SIRts(ts) = 29.5344;
        numts(3) = numts(3) + 1;
    elseif SIRdb(ts) > -2.04 && SIRdb(ts) < .66
        SIRts(ts) = 47.502;
        numts(4) = numts(4) + 1;
    elseif SIRdb(ts) > .66 && SIRdb(ts) < 2.84
        SIRts(ts) = 75.8016;
        numts(5) = numts(5) + 1;
    elseif SIRdb(ts) > 2.84 && SIRdb(ts) < 4.73
        SIRts(ts) = 110.502;
        numts(6) = numts(6) + 1;
    elseif SIRdb(ts) > 4.73 && SIRdb(ts) < 6.38
        SIRts(ts) = 148.1508;
        numts(7) = numts(7) + 1;
    elseif SIRdb(ts) > 6.38 && SIRdb(ts) < 8.78
        SIRts(ts) = 186.0516;
        numts(8) = numts(8) + 1;
    elseif SIRdb(ts) > 8.78 && SIRdb(ts) < 11.49
        SIRts(ts) = 241.1766;
        numts(9) = numts(9) + 1;
    elseif SIRdb(ts) > 11.49 && SIRdb(ts) < 13.27
        SIRts(ts) = 303.1938;
        numts(10) = numts(10) + 1;
    elseif SIRdb(ts) > 13.27 && SIRdb(ts) < 16.52
        SIRts(ts) = 344.043;
        numts(11) = numts(11) + 1;
    elseif SIRdb(ts) > 16.52 && SIRdb(ts) < 19.71
        SIRts(ts) = 418.6098;
        numts(12) = numts(12) + 1;
    elseif SIRdb(ts) > 19.71 && SIRdb(ts) < 23.12
        SIRts(ts) = 491.6898;
        numts(13) = numts(13) + 1;
    elseif SIRdb(ts) > 23.12 && SIRdb(ts) < 26.37
        SIRts(ts) = 569.9484;
        numts(14) = numts(14) + 1;
    elseif SIRdb(ts) > 26.37 && SIRdb(ts) < 28.79
        SIRts(ts) = 644.5152;
        numts(15) = numts(15) + 1;
    elseif SIRdb(ts) > 28.79
        SIRts(ts) = 699.8922;
        numts(16) = numts(16) + 1;
    end
end

propnumts = numts./tamcenu;
xts = 0:15;

figure(5)
%     hist(SIRts,15);
    bar(xts,propnumts)
    title('Distribución de probabilidades de las tasas obtenidas, k=1');
    grid on;
    ax = gca;
    ax.XTick = 0:15;
    ax.XTickLabels = {'0','19.1898','29.5344','47.5020','75.8016','110.5020','148.1508','186.0516','241.1766','303.1938','344.0430','418.6098','491.6898','519.9484','644.5152','699.8922'};
    ax.XTickLabelRotation = 90;
    xlabel('Tasa asignada');
    ylabel('Probabilidad')

%Parte2--------------------------------------------------------------------
Mx = [1,2,10,20,50];
% for im = 1:length(Mx)
    im=1;
    usuarios = 50;
    varTiempo = 1; %Incrementar 1 en 1 suponiendo que es 1 ms
    rjtast = SIRts(1:usuarios); %rj*(t)
    M = Mx(3);
    Tval = 1000;
    Delta = 1;
    % Rjtint = zeros(1,usuarios); 
    Rjt = zeros(1,usuarios); 
    rjt = zeros(1,usuarios); 
    rsum = zeros(1,usuarios); 
    tasapromedio = zeros(1,usuarios);
    for T = 0:Tval-1
        if(T<M)
            Rjtpasada = Rjt;
            if T == 0
                Rjt = rsum./(0.0000001);
            else
                Rjt = rsum./(T);
            end
        else
            Rjtpasada = Rjt;
            Rjt = ( ((M-1)/(M)).*Rjtpasada ) + ( (1/M).*rsumpasada );
        end

        for t = 1:usuarios
            kj(t) =  rjtast(t)./(Delta + Rjt(t));
        end

        [m n] = max(kj);
        rjt = zeros(1,usuarios);
        rjt(n) = m; 
        tasapromedio(n) = tasapromedio(n) + m;

        ri(T+1) = rjt(n);
        rsumpasada = rsum;
        rsum = rsum + rjt;
    end
    tjpromedio = tasapromedio./Tval;
    tjpromedio2 = tjpromedio.^2;
    J = ((sum(tjpromedio)).^2)./(tjpromedio2);

    figure(im+5)
    subplot(1,2,1);
    plot(1:M,ri(1:M))
    title(' Grafica ')
    xlabel('M');
    ylabel('r(kb/s)');

    subplot(1,2,2);
    plot(1:M,J(1:M))
    title(' Grafica ')
    xlabel('M');
    ylabel('J');
% end

% Rjt(t) = ( ((M-1)/(M))*Rjt(t-1) ) + ( (1/M)* SIRast_jt(t) )

