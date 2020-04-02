clear all;
clc;
%1.- Simulacion de la transmisión ----------------------------------------- 
t = 0:1:15;
posw2 = [-1 -1;
         -1 +1];  
posw4 = [posw2 posw2; 
         posw2 -posw2];     
posw8 = [posw4 posw4;
         posw4 -posw4];
posw16= [posw8 posw8; 
         posw8 -posw8];
cxyt = posw16;

%Para generar cada usuario
% for i = 1:16
%     numal = rand;
%     if numal < 0.5
%         x(1,i) = -1;
%     elseif numal > 0.5  
%         x(1,i) = 1;
%     end
% end
x=[1 -1 -1 1 -1 1 1 -1 -1 -1 1 -1 -1 -1 -1 -1];

for i = 1:16
    xc(i,:) = cxyt(i,:).*x(1,i);
end

s1d = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
for i = 1:16
        s1d(1,:) = s1d(1,:) + xc(i,:);
end

%2.- Simulacion de la recepción ------------------------------------------- 
for i = 1:16
    multscxyt(i,:) = s1d(1,:) .* cxyt(i,:);
end

s2b = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
trans = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
for i = 1:16 %filas
    for ii = 1:16
            s2b(i) = s2b(i) + multscxyt(i,ii);
    end
    if  s2b(i) < 0 
        trans(i) = -1;
    elseif s2b(i) > 0 
        trans(i) = 1;
    end
end

figure(1)
subplot(3,1,1);
stem(t,x);
grid on;
xlabel('t');
ylabel('amplitud');
title(' ');

subplot(3,1,2);
stem(t,sgauss);
grid on;
xlabel('t');
ylabel('amplitud');
title(' ');

subplot(3,1,3);
stem(t,trans);
grid on;
xlabel('t');
ylabel('amplitud');
title(' ');

%Trans
trans
%Usuario
x 

% % %3.- Simulacion con ruido blanco aditivo y gaussiano----------------------- 
% for i = 1:16
%     r(i) = random('Normal',0,2);
%     sgauss(1,i) =  s1d(1,i) + r(i);  
% end
% 
% s1d
% sgauss
