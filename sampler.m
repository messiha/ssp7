clear;clc;
%Source Motion Model Parameters
delta_t = 0.375;
N_samples = 1000;
N_steps = int8(10/delta_t);
beta = 2;
v_bar = 1;
a = exp(-beta*delta_t);
b = v_bar*sqrt(1-a^2);
F = [eye(2), a*delta_t*eye(2);zeros(2,2),a*eye(2)];
Q = [b^2*delta_t^2*eye(2),zeros(2,2);zeros(2,2),b^2*eye(2)];



%True source and microphone motion models
source = zeros(N_steps+1,4);
source(1,1:2) = 0.5*[randi([3 9],1,2)] ;
mic = zeros(N_steps+1,2,3);
mic(1,1,:) = [1.35,1,1.5];
mic(1,2,:) = [1.65,1,1.5];
mic(end,1,:)  = [4.85,1,1.5];
mic(end,2,:) = [5.15,1,1.5];

for i = 1:3
	mic(:,1,i) = linspace(mic(1,1,i),mic(end,1,i),N_steps+1);
	mic(:,2,i) = linspace(mic(1,2,i),mic(end,2,i),N_steps+1);
end

%rir parameters and model
room_dim = [6 6 2.5];
fs = 8000;
c = 340;
rt60  = 0.5;
rir_samples = fs*delta_t;
rir = zeros(N_steps+1,2,rir_samples);
%rir(1,:,:) = rir_generator(c,fs,mic(1,:,:),[source(1,1:2) 0],room_dim,rt60,rir_samples);

%Sampled Source Positions
source_samp = zeros(N_steps+1,N_samples,4);
source_samp(1,:,:) = [0.5*randi([3 9],N_samples,2) zeros(N_samples,2)];


%Grid points
[X,Y] = meshgrid(0:0.1:6,0:0.1:6);
grid_pts = [X(:),Y(:)];
kdeprob = zeros(N_steps,size(grid_pts,1));

%Signal Array : sig 

%STFT params
window_size = 100;
window = rectwin(window_size);



for t = 1:N_steps
	temp_source = mvnrnd(F*reshape(source(t,:),[4 1]),Q,1);
    a = temp_source(1:2) <0 ;
    b = temp_source(1:2) > 6;
	while ((sum(a(:)) > 0) || (sum(b(:)) > 0))
			temp_source = mvnrnd(F*reshape(source(t,:),[4 1])   ,Q,1);
    end
    
	source(t+1,:) = temp_source;		
	rir(t,:,:) = rir_generator(c,fs,reshape(mic(t,:,:),[2,3]),[source(t,1:2) 0],room_dim,rt60,rir_samples);	
	
	%Convolve with signal from TIMIT database
	%S1 = conv(rir(t,:,1),signal(t,:))
	%S2 = conv(rir(t,:,2),signal(t,:))
	%STFT with signal array
	%Z(t,1,:) = spectrogram(S1,window,[],[],fs)
	%Z(t,2,:) = spectrogram(S2,window,[],[],fs)

	for j = 1:N_samples
		
		temp_samp = mvnrnd(F*reshape(source_samp(t,j,:),[4,1]),Q,1);
        a = temp_samp(:,1:2) < 0;
        b =  temp_samp(:,1:2) > 6;
		while ( (sum(a(:)) > 0) || (sum(b(:)) > 0))
			temp_samp = mvnrnd(F*reshape(source_samp(t,j,:),[4,1]),Q,1);
        end
   		source_samp(t+1,j,:) = temp_samp;
		%Update equations for weight
		%w(t+1,j,:) = w(t,j,:)
    end
    
	%kdeprob(i,:),pts = ksdensity(reshape(source_samp(t,:,1:2),[N_samples,2]),grid_pts,'Weights',reshape(w(t,:,:),[N_samples,2]));
end    
contour(X,Y,reshape(kdeprob(1,:),size(X)));
		
