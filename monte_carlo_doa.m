clear all;
%---------------------SETUP MODEL------------------------------------------
%uniform linear array of 11 sensors with isotropic power kept at 0.5 lambda

c = physconst('LightSpeed');
fc = physconst('LightSpeed');              % Operating frequency
lambda = c/fc;  %wavelength of 1 m
uniform_distance = lambda/2;

ula = phased.ULA('NumElements',11,'ElementSpacing',uniform_distance);

%---------------------signal model-----------------------------------------
ang1 = [-20; 0];          % First signal
ang2 = [30; 0];          % Second signal


%signal model with 10 sources at equally spaced between -60 and +60
angs = [ang1 ang2];
ang = [ang1(1,1) ang2(1,1)]/180*pi;

pos = getElementPosition(ula)/lambda; 

Nsamp = 100;        % 100 samples
nPower = 1/10;       % SNR of 10dB

% Assignment Matrix
A = zeros(2,11);
for k = 1:2
    A(k,:) = exp(-1i*2*pi*lambda/2*sin(ang(k))/lambda*[0:10]);
end    
A = transpose(A);
    
rs = rng(1996);     % pseudo random states for reproducible results
signal = sensorsig(pos,Nsamp,angs,nPower); %generated signal model samples

%-----------------------BEAMSCAN-------------------------------------------

spatialspectrum = phased.BeamscanEstimator('SensorArray',ula,...
            'OperatingFrequency',fc,'ScanAngles',-90:90);
spatialspectrum.DOAOutputPort = true;
spatialspectrum.NumSignals = 2;
[~, ang_beamscan] = step(spatialspectrum, signal)
% plotSpectrum(spatialspectrum);
% hold on;
%--------------------MVDR--------------------------------------------------

mvdrspatialspect = phased.MVDREstimator('SensorArray',ula,...
        'OperatingFrequency',fc,'ScanAngles',-90:90,...
        'DOAOutputPort',true,'NumSignals',2);
[~,ang_mvdr] = step(mvdrspatialspect, signal)
% plotSpectrum(mvdrspatialspect);

%------------------MUSIC---------------------------------------------------

musicspatialspect = phased.MUSICEstimator('SensorArray',ula,...
        'OperatingFrequency',fc,'ScanAngles',-90:90,...
        'DOAOutputPort',true,'NumSignalsSource','Property','NumSignals',2);
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
    yig(iter,1) = 1/pencil;
end

[val, ind] = sort(yig, 'descend');
ang_ig = ind(1:2)'-91



plot(phi, yig);
plot(ang/pi*180, zeros(length(ang)), 'k*')
legend('MVDR','MUSIC','IG Pencil', 'TRUE DOA');



    
    
        
        
    







