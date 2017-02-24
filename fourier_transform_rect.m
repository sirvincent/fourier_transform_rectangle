%% Fourier transform of a rectangle function in 2D and 3D
% showcasing the usage of the discrete fourier transform and the extra steps
% needed to make it coincide with the analytical result of fourier
% transforming a rectangle. 
% For more information about fourier transforms and the 
% discrete fourier transforms see:

% https://en.wikipedia.org/wiki/Fourier_transform
% https://en.wikipedia.org/wiki/Discrete_Fourier_transform

% created by Vincent Post (2016)
%% Close all graphs and clear the memory. 
close all;
clear all;

%% Input of user defined parameters
% N = number of samples
% L = support width in real space
% b = width of rectangle function.
N = 128;
L = 12;
b = 1;

%% Creation rectangle and fourier transfom of rectangle
deltax = L/N;   % sample spacing 
x = -L/2:deltax:L/2-deltax; % spatial coordinates
func = double(abs(x)<= b/2); % rectangle

% Setting up frequency domain values (with an even number of elements)
deltaFreq = 1/L;
freq = -1/(2*deltax):deltaFreq:1/(2*deltax)-deltaFreq;
% apply fft
funcfourier = fft(func);

% the analytical result of the fourier transform of a rectangle is a 
% sinc function (see also Goodman: fourier optics)
% Create analytical sinc function for comparison
func1Dsinc = b*sinc(b*freq);

%% Plotting of rectangle and fft rect which is compared to the sinc
%       by comparing absolute, imaginary and real part seperately
subplot(2,2,1);
plot(x,func,'-or');
axis([-L/2, L/2, 0.0, 1.2]);
xlabel('position');
ylabel('function value');
title('rectangle');

% Compare absolute value fft rect with the analytical sinc
subplot(2,2,2);
plot(freq,abs(funcfourier),'-or'); 
hold on;
plot(freq,abs(func1Dsinc),'-ob');
xlabel('frequency');
ylabel('function value');
legend('fft rect', 'sinc');
title('absolute value');

% Compare real value fft rect with analytical sinc
subplot(2,2,3);
plot(freq,real(funcfourier),'-or');
hold on;
plot(freq,real(func1Dsinc),'-ob');
xlabel('frequency');
ylabel('function value');
legend('fft rect', 'sinc');
title('real value');

% Compare imaginary value fft rect with analytical sinc
subplot(2,2,4);
plot(freq,imag(funcfourier),'-or');
hold on;
plot(freq,imag(func1Dsinc),'-ob');
xlabel('frequency');
ylabel('function value');
legend('fft rect', 'sinc');
title('imaginary value');


%% Interpertation fft plots
%{ 
    The absolute value of the fft of rect is in comparison to the absolute value of
    sinc split and flipped in the middle. Also, the amplitude is different.
    The same holds for the real value of the fft of the rect but with an extra
    sinusoidal factor...the sinusoidal factor may be attributed to the fact that
    seen by the fft the array is shifted. Hence the fft obtains a phase factor
    according to the shift theorem (See wiki or the book of Goodman Fourier Optics).
    The imaginary value of the fft of the rect is the same as the imaginary value of a
    sinc, just a 0 line. I believe this is due that the rect
    has a fixed phase.
%}

%% Shifting the fourier transform
% To solve this problem that the fft of the rectangle is flipped. Due to the
% implementation of the discrete fourier transform. A function called
% fftshift is implemented in matlab to shift the DC component (0 frequency
% component) to the middle of the array
funcfourier = fftshift(funcfourier);    
%% Comparing rectangle and sinc after fftshift
% absolute part
figure;
subplot(2,2,1);
plot(freq,abs(funcfourier),'-or');
hold on;
plot(freq,abs(func1Dsinc),'-ob');
xlabel('frequency');
ylabel('function value');
legend('fft rect', 'sinc');
title('absolute value');

% real part
subplot(2,2,2);
plot(freq,real(funcfourier),'-or');
hold on;
plot(freq,real(func1Dsinc),'-ob');
xlabel('frequency');
ylabel('function value');
legend('fft rect', 'sinc');
title('real value');

% imaginary value
subplot(2,2,3);
plot(freq,imag(funcfourier),'-or');
hold on;
plot(freq,imag(func1Dsinc),'-ob');
xlabel('frequency');
ylabel('function value');
legend('fft rect', 'sinc');
title('imaginary');

%% Interpertation fftshift plots
%{ 
    The real and absolute graphs are now flipped, creating a peak in the
    middle, looking more like the original sinc.
    Even though, the amplitude is still off and for the real part there are
    still fluctuations.
    The imaginary part remains the same as before, as good as 0
%}

%% Correcting the fourier transform
% Due to the discretization of the discrete fourier transform we need to
% multiply the fourier transform with the sampple spacing to let the peak
% positions and amplitude coincide. Also the sinusoidal term in the real part, 
% arising from the fact that the fft seems the array shifted, are removed.
% This is done by introducing an extra phase factor in the exponential that cancels the
% obtained phase factor from the shift theorem
corvec = deltax*exp(-1i*2*pi*freq*L/2);
funcfourier = corvec.*funcfourier;

%% Comparing results with the correction vector
% absolute part
figure;
subplot(2,2,1);
plot(freq,abs(funcfourier),'-or');
hold on;
plot(freq,abs(func1Dsinc),'-ob');
xlabel('frequency');
ylabel('function value');
legend('fft rect', 'sinc');
title('absolute value');

% real part
subplot(2,2,2);
plot(freq,real(funcfourier),'-or');
hold on;
plot(freq,real(func1Dsinc),'-ob');
xlabel('frequency');
ylabel('function value');
legend('fft rect', 'sinc');
title('real value');

% imaginary value
subplot(2,2,3);
plot(freq,imag(funcfourier),'-or');
hold on;
plot(freq,imag(func1Dsinc),'-ob');
xlabel('frequency');
ylabel('function value');
legend('fft rect', 'sinc');
title('imaginary');

%% interpertation corvec plots
%{ 
    The fft is almost exactly the same as the analytical solution.
    The lines do not precisely coincide which might be because of the sampling of
    the fourier transform.
%}

%% inverse fourier transforming(ifft) to obtain the rectangle back
invfuncfourier = ifftshift(ifft(funcfourier))./corvec;

%% Plotting inverse fourier transform
% absolute value
figure;
subplot(2,2,1);
plot(x,abs(invfuncfourier),'-or');
hold on;
plot(x,abs(func),'-ob');
xlabel('x');
ylabel('function value');
legend('ifft rect', 'original rect');
title('absolute value');

% real value
subplot(2,2,2);
plot(x,real(invfuncfourier),'-or');
hold on;
plot(x,real(func),'-ob');
xlabel('x');
ylabel('function value');
legend('ifft rect', 'ooriginal rect');
title('real value');

% imag value
subplot(2,2,3);
plot(x,imag(invfuncfourier),'-or');
hold on;
plot(x,imag(func),'-ob');
xlabel('x');
ylabel('function value');
legend('ifft rect', 'original rect');
title('imaginary value');
%}
%% Interpertation inverse fourier transform (ifft) plots
%{ 
    The corvec applied in the frequency domain needs to be removed again.
    Also the array needs to be shifted again, now by using ifftshift.
    As a consequence we see the original result coinciding with the ifft
%}

%% fourier transforming in 3D of the rectangle function
% setting up rectangle and real space variables.
deltax = L/N;
x = -L/2:deltax:L/2-deltax;
y = -L/2:deltax:L/2-deltax;
[XX,YY] = meshgrid(x,y);
func2D = double(abs(XX)<= b/2).*double(abs(YY)<= b/2);

% fourier transform and setting up the corresponding frequency variables
% including the correction vector.
deltaFreq = 1/L;
freqx = -1/(2*deltax):deltaFreq:1/(2*deltax)-deltaFreq;
freqy = -1/(2*deltax):deltaFreq:1/(2*deltax)-deltaFreq;
[FREQX, FREQY] = meshgrid(freqx, freqy);
corvec = deltax^2*exp(-1i*2*pi*(FREQX+FREQY)*L/2);
funcfourier2D = corvec.*fftshift(fft2(func2D)); % using matlabs fft2 not fft

func2Dsinc = sinc(FREQX).*sinc(FREQY);  % set b = 1

%% Plot 3D fourier transform of rectangle
% rectangle
figure;
subplot(2,2,1);
surfc(XX,YY,func2D);
colormap default;
view([-37,36]);
xlabel('x');
ylabel('y');
zlabel('function value');
title('rect'); 

% absolute value
subplot(2,2,2);
surfc(FREQX,FREQY, abs(funcfourier2D));
hold on;
surfc(FREQX,FREQY, abs(func2Dsinc)); 
colormap default;
view([-37,36]);
xlabel('freq x');
ylabel('freq y');
zlabel('function value');
title('absolute value fft rect'); 

% real value
subplot(2,2,3);
surfc(FREQX,FREQY, real(funcfourier2D));
%surfc(FREQX,FREQY, real(func2Dsinc)); 
colormap default;
view([-37,36]);
xlabel('freq x');
ylabel('freq y');
zlabel('function value');
title('real value fft rect'); 

% imaginary value
subplot(2,2,4);
surfc(FREQX,FREQY, imag(funcfourier2D)); 
%surfc(FREQX,FREQY, imag(func2Dsinc)); 
colormap default;
view([-37,36]);
xlabel('freq x');
ylabel('freq y');
zlabel('function value');
title('imaginary value fft rect'); 


%% Interpertation 3D fourier transform plots
%{ 
    The fft of rectangle in 3D looks as a sinc in 3D(not shown when script is 
    run but it is calculated).
    Beautifull! :)
%}