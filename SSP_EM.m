% almost completed
function prob = SSP_EM(mic_pos,particle_pos,z,epsilon)

%%
% t : a particular time instant we are working on
% mic_pos =2X3
% particle_pos = 2 x J

prob = ones(1,J);
c = 3*10^8;

% z = 2xK
Ts = 1/8000;
K = size(z);
K = K[2];                    %K : Frequency bins count

J = size(particle_pos);
J = J[2];                   % J : Sampled positions for source

    
% h and tau are given
tau = zeros(K,2,2);
for k = 1:K
    xx = sinc(2*pi*k*norm(mic_pos(1,:)-mic_pos(2,:))/(K*Ts*c));
    for ii = 1:2
         for jj = 1:2
             tau(k,ii,jj) = xx;
             if ii==jj
                tau(k,ii,jj) = xx + epsilon;
             end
         end
    end
end


h = zeros(J,K,2);
for jj = 1:J
    for k = 1:K
        for ii = 1:2
            c = mic_pos(1,:)/2+mic_pos(2,:)/2;
            %a = atan()   % to ask
            h(jj,k,ii) = exp(2j*pi*k*norm(mic_pos(ii,:)-c)*cos(gamma)/(K*Ts*c));   % assuming centre is the reference microphone then gamma_m = 0 or 180
        end
    end
end

phi_y = zeros(J,K,2,2);
phi = zeros(J,K,2,2);     % _init_ phi, sai, phi_r and phi_y as zeros tensors
phi_r = zeros(J,K,2,2);   % phi, phi_y and phi_r are 2X2 matrices; phi(t,k) is a 2X2
sai = ones(J,1)/J;        % for each j and k

shift = 10000;  % a big number
iter = 0;       % counter
epsilon = 0.0001; % percision
                   
%formatSpec = 'iteration: %d, error: %2.4f, mu1: [%2.4f %2.4f], mu2: [%2.4f %2.4f] \n';
%%
while shift > epsilon
    iter = iter + 1;
    if iter == 1
        mu = ones(J,K)/J;              % mu : is an J X K dimen matrix
    else
        mu1 = ones(J,K)/J;
        for k = 1:K
            for jj = 1:J
                p = mvnpdf(z(:,k),[0;0],reshape(phi(jj,k,:,:),[2,2]));
                mu1(jj,k) = sai(jj,1)*p;    % completed it with reshape(phi(jj,k,:,:),[2,2])
            end
            mu1(:,k) = mu(:,k)/sum(mu(:,k));
        end
    end
    
    phi_y1 = zeros(J,K,2,2);
    phi1 = zeros(J,K,2,2);     % _init_ phi, sai, phi_r and phi_y as zeros tensors
    phi_r1 = zeros(J,K,2,2);   % phi, phi_y and phi_r are 2X2 matrices; phi(t,k) is a 2X2
    sai1 = ones(J,1)/J;        % for each j and k
    
    for k = 1:K
        for jj = 1:J
            sai1(jj) = mean(mu1(jj,:));
            b(jj,:) = inv(reshape(tau(k,:,:),[2,2]))*reshape(h(jj,k,:),[2,1])/(reshape(h(jj,k,:),[2,1])'*inv(reshape(tau(k,:,:),[2,2]))*reshape(h(jj,k,:),[2,1]));
            phi_r1(jj, k,:,:) = (z(:,k))*(eye(2) - b(jj,:)'*(h(jj,k,:))')*inv(reshape(tau(k,:,:),[2,2]))*(z(:,k))';
            phi_y1(jj, k,:,:) = b(jj,:)*((z(:,k))*z(:,k)'-reshape(phi_r1(jj,k,:,:),[2,2])*reshape(tau(k,:,:),[2,2]))*(b(jj,:))';
            phi1(jj,k,:,:) = reshape(h(jj,k,:),[2,1])*reshape(h(jj,k,:),[2,1])'*reshape(phi_y1(jj,k,:,:),[2,2]) + reshape(tau(k,:,:),[2,2])*reshape(phi_r1(jj,k,:,:),[2,2]);  
        end
              
    end
    
    
    if iter>1
        shift = calc_distance(phi,phi1,mu,mu1,J,K);
        phi = phi1;
        mu = mu1;
    end
end

for jj = 1:J                %confirm
    disp('1M US $')
%     for k = 1:K
%         prob(1,jj) = prob(1,jj)*sai
%     end
end

end

