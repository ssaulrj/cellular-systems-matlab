clear all;
clc;
close all;
%--------------------------------------------------------------------------
Tu = 0.001;
A = 1e7;
Ts = 1/A;
t = [0:(Ts):Tu - Ts];
Df = 1000;

% Para generar cada usuario
for i = 1:16
    numal = rand;
    if numal < 0.5
        ax(1,i) = -1;
    elseif numal > 0.5  
        ax(1,i) = 1;
    end
end

ax1s = ones(1,16);

Nc = 16;
for kk = 1:Nc  
    if ax1s(1,kk) == 1
        a2 = ones(1,10000);
    end
    if ax(1,kk) == 1
        a = ones(1,10000);
    elseif ax(1,kk) == -1
        a = -ones(1,10000);
    end  
    c = exp(j.*2.*pi.*(kk-1).*Df.*t);
    c2 = exp(j.*2.*pi.*(kk-1).*Df.*t);
    sm = a(1,kk).*c;
    sm2 = a2(1,kk).*c2;
    smx(kk,:) = sm;
    smx2(kk,:) = sm2;
end

%Alternativa es usar la funcion SUM(X)
sumax = zeros(1,10000);
sumax2 = zeros(1,10000);
for i= 1:Nc
    sumax(1,:) = sumax(1,:)+smx(i,:);
    sumax2(1,:) = sumax2(1,:)+smx2(i,:);
end

% figure(1)
% for i = 1:16 
%     subplot(4,4,i)
%     plot(real(smx(i,:)));
%     xlabel('t');
%     ylabel('Amplitud');
%     grid on;
% end

N = 256;
xespectro = N .* ifft(ax,N);
[m n] = size(xespectro);
nn = [0:n-1];
tn = [0:Tu/N:Tu-(Tu/N)];

xespectro2 = N .* ifft(ax1s,N);
[m2 n2] = size(xespectro2);
nn2 = [0:n2-1];
tn2 = [0:Tu/N:Tu-(Tu/N)];

xt = real(sumax);
xtfft = fftshift(fft(xt, N));

xt2 = real(sumax2);
xtfft2 = fftshift(fft(xt2, N));

% nfft = 1024;
% y = fft(xt,nfft);
% y = y(1:nfft/2);
% mi = abs(y).^2;
% % f = (0:length(y)-1)*A/length(y);
% f = (0:nfft/2-1)*A/nfft;
% 
% % xt2=real(sumax2);
%
% y2 = fft(xt2,nfft);
% y2 = y2(1:nfft/2);
% mi2 = abs(y2).^2;
% % f = (0:length(y)-1)*A/length(y);
% f2=(0:nfft/2-1)*A/nfft;

figure(2)
subplot(2,1,1);
plot(t,xt);
xlabel('t');
ylabel('Amplitud');
title('Señal OFDM');
grid on;
hold on;
stem(tn,real(xespectro));
xlabel('');
ylabel('Amplitud');
title('Espectro ifft');
grid on;

subplot(2,1,2);
hold on;
plot(abs(xtfft));
xlabel('f');
ylabel('Amplitud');
title('fft x(t)');
% xlim([-.5 1.5]);
grid on;

figure(3)
subplot(2,1,1);
grid on;
plot(t,xt2);
xlabel('t');
ylabel('Amplitud');
title('Señal OFDM a_k=unos');
hold on;
stem(tn2,real(xespectro2));
xlabel('t');
ylabel('Amplitud');
title('Espectro ifft a_k=unos');


subplot(2,1,2);
hold on;
grid on;
plot(abs(xtfft2));
xlabel('f');
ylabel('Amplitud');
title('fft x(t)');

