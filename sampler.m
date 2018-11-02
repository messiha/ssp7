delta_t = 0.375;
N_samples = 1000;
N_steps = int8(10/delta_t);
%Source motion Model
beta = 2;
v_bar = 1;
a = exp(-beta*delta_t);
b = v_bar*sqrt(1-a^2);
F = [eye(2), a*delta_t*eye(2);zeros(2,2),a*eye(2)];
Q = [b^2*delta_t^2*eye(2),zeros(2,2);zeros(2,2),b^2*eye(2)];

%Initial Samples
S = zeros(N_steps+1,N_samples,4);
S(1,:,:) = 0.5*randi([3 9],N_samples,4);
[X,Y] = meshgrid(0:0.1:6,0:0.1:6);
pts = [X(:),Y(:)];
kdeprob = zeros(N_steps,size(pts)(1))
for t = 1:N_steps
	for j = 1:N_samples
		S(t+1,j,:) = mvnrnd(F*reshape(S(t,j,:),[N_samples,4]),Q,1);
		%Update equations for weight
		%w(t+1,j,:) = w(t,j,:)
	kdeprob(i,:),pts = ksdensity(reshape(S(t,:,1:2),[N_samples,2]),pts,'Weights',reshape(w(t,:,:),[N_samples,2]));

contour(X,Y,reshape(kdeprob(1,:),size(X));
		
