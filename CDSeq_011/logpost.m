function y = logpost(WD,T,theta,phi,Bet,Alp)
% coder : Kai Kang

[tr,tc]=size(theta);
[pr,pc]=size(phi);
if T~=tr
    error('Error: theta should be a T by samplesize matrix');
end

if T~=pc 
    error('Error: phi should be a gene by T matrix');
end
mu=1e8;
t1 = sum(log(mu*phi),2)'*(Bet-1)';
t2 = sum(log(theta),2)'*(Alp-1)';
t3 = sum(sum(WD'.*log(mu*theta'*phi')));
y = t1 + t2 + t3;

end