clear all;
%---------------------SETUP MODEL------------------------------------------
%uniform linear array of 11 sensors with isotropic power kept at 0.5 lambda

c = physconst('LightSpeed');
fc = physconst('LightSpeed');              % Operating frequency
lambda = c/fc;  %wavelength of 1 m
uniform_distance = lambda/2;

ula = phased.ULA('NumElements',11,'ElementSpacing',uniform_distance);

%---------------------signal model-----------------------------------------
ang1 = [-60; 0];          % First signal
ang2 = [-47; 0];          % Second signal
ang3 = [-34; 0];
ang4 = [-20; 0];
ang5 = [-7; 0];
ang6 = [6; 0];
ang7 = [20; 0];
ang8 = [33; 0];
ang9 = [46; 0];
ang10 = [60; 0];

%signal model with 10 sources at equally spaced between -60 and +60
angs = [ang1 ang2 ang3 ang4 ang5 ang6 ang7 ang8 ang9 ang10 ];
ang = [ang1(1,1) ang2(1,1) ang3(1,1) ang4(1,1) ang5(1,1) ang6(1,1)...
    ang7(1,1) ang8(1,1) ang9(1,1) ang10(1,1)]/180*pi;

pos = getElementPosition(ula)/lambda; 

Nsamp = 100;        % 100 samples
nPower = 1/10;       % SNR of 10dB

% Assignment Matrix
A = zeros(10,11);
for k = 1:10
    A(k,:) = exp(-1i*2*pi*lambda/2*sin(ang(k))/lambda*[0:10]);
end    
A = transpose(A);
    
rs = rng(1996);     % pseudo random states for reproducible results
signal = sensorsig(pos,Nsamp,angs,nPower); %generated signal model samples

%-----------------------BEAMSCAN-------------------------------------------

spatialspectrum = phased.BeamscanEstimator('SensorArray',ula,...
            'OperatingFrequency',fc,'ScanAngles',-90:90);
spatialspectrum.DOAOutputPort = true;
spatialspectrum.NumSignals = 10;
[~, ang_beamscan] = step(spatialspectrum, signal)
% plotSpectrum(spatialspectrum);
% hold on;
%--------------------MVDR--------------------------------------------------

mvdrspatialspect = phased.MVDREstimator('SensorArray',ula,...
        'OperatingFrequency',fc,'ScanAngles',-90:90,...
        'DOAOutputPort',true,'NumSignals',10);
[~,ang_mvdr] = step(mvdrspatialspect, signal)
% plotSpectrum(mvdrspatialspect);

%------------------MUSIC---------------------------------------------------

musicspatialspect = phased.MUSICEstimator('SensorArray',ula,...
        'OperatingFrequency',fc,'ScanAngles',-90:90,...
        'DOAOutputPort',true,'NumSignalsSource','Property','NumSignals',10);
[~,ang_music] = musicspatialspect(signal)
ymvdr = mvdrspatialspect(signal);
ymusic = musicspatialspect(signal);
helperPlotDOASpectra(mvdrspatialspect.ScanAngles,...
        musicspatialspect.ScanAngles,ymvdr,ymusic,'ULA')
hold all;

%----------------------IG Pencil-------------------------------------------

% signal = transpose(signal);

%Since true covariance matrix is unavailable it is replaced by sample
%covariance matrix
R_hat = signal'*signal;

%eigen values and eigen vectors of sample covariance matrix
[N,V] = eig(R_hat);

%Linear Search for maximas from -pi/2 to pi/2
phi = -90:1:90;
yig = zeros(length(phi),1);
for iter=1:length(phi)
    a_phi = zeros(11,1);
    for j = 1:11
        a_phi(j,1) = exp(-1i*2*pi*(j-1)*uniform_distance/lambda*...
            sin(phi(1,iter)/180*pi));
    end
    pencil = abs(a_phi'*inv(R_hat)*a_phi);
    yig(iter,1) = pencil;
end
yig = yig.^-1;
inverse_func = yig;
yig = yig/max(yig);
yig =5*log10(yig);
% [val, ind] = sort(yig, 'descend');
% ang_ig = ind(1:10)-91

plot(phi, yig);
plot(ang/pi*180, zeros(length(ang)), 'k*')
legend('MVDR','MUSIC','IG Pencil', 'TRUE DOA');



    
    
        
        
    







