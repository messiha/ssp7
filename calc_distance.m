%  completed


function d = calc_distance(phi,phi1,mu,mu1,J,K)


d = max(max(abs(mu-mu1))');

for jj = 1:J
    for k = 1:K
        dd = max(max(abs(reshape(phi(jj,k,:,:),[2,2])-reshape(phi1(jj,k,:,:),[2,2])))');
        if dd > d
            d = dd;
        end
    end
end

end