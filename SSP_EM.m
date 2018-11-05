% J : Sampled positions for source   %K : Frequency bins count %t : a particular time instant we are working on
% phi_y and phi_r are 2X2 matrices; phi(t,k) is a 2X2   %L : no of iterations performed for EM

% mu : is an J X K dimen matrix
% sai : J dim vector

% initialize sai , maybe phi_r and phi_y stored in struct params
% assume h and tau are given

% params.sai = 


shift = 10000;  % a big number
iter = 0;       % counter
epsilon = 0.001; % percision
formatSpec = 'iteration: %d, error: %2.4f, mu1: [%2.4f %2.4f], mu2: [%2.4f %2.4f] \n';

while shift > epsilon
    iter = iter + 1;
    mu = zeros(J,K);
 
    for k = 1:K
        for jj = 1:J 
            phi(jj,k) = h(jj,k,:)*(h(jj,k,:))'phi_y(jj, k,:) + tau(k,:,:)*phi_r(jj, k);
            mu(jj,k) = sai(jj)*prob_nor(z(k));
        end
        mu(:,k) = mu(:,k)/sum(mu(:,k));
    end
    
    for k = 1:K
        for jj = 1:J
            sai(jj) = mean(mu(jj,:));
            b(jj,:) = inv(tau(k,:,:))*h(jj,k,:)/(h(jj,k,:)'inv(tau(k,:,:))*h(jj,k,:);
        end
        
        for jj = 1:J            
            phi_r(jj, k) = (z(k,:))*(eye(2) - b(jj,:)'*(h(jj,k,:))')*inv(tau(k,:,:))*(z(k,:))';
            phi_y(jj, k) = b(jj,:)*((z(k,:))'*z(k,:)-phi_r(jj, k)*tau(k,:,:))*(b(jj,:))';
        end
    end
      
    
    shift = calc_distance(Param, Param_);
end

