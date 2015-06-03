function kp = Kmore_Prob(pV,k,flag)

%compute k-more prob: at least k events will occur;
%2008-4-27, refine;
%fix the bug: thanks to Kensuke Onuma
if nargin<3
    flag = 1;
end
if nargin<2
    k = 0;
end
[m,n] = size(pV);

if flag==1
    if k==n
        kp = Com_AND(pV);
    elseif k==1
        kp = Com_OR(pV);
    else
        for i=1:n
            %%should be a better way,- do not need to store all A{i}, only
            %%two matrices should be enough
            A{i} = zeros(m,n);
            A{i}(:,1) = Com_OR(pV(:,[1:i]));
            for k0=2:min(i,k)                
                A{i}(:,k0) = A{i-1}(:,k0-1).*pV(:,i)+A{i-1}(:,k0).*(1-pV(:,i));                
            end
        end
        kp = A{n}(:,k);   
    end
else%by order statistic
    tmp = sort(pV,2);
    kp = tmp(:,end-k+1);
end



% function nP = Addone(tP,p)
% 
% %tP is 1xlen, tP(i) denotes the prob at least (i-1) events will occue
% %p is the prob of a new event
% 
% nP = [0 p*tP]+[tP 0];


function kp = Com_AND(pV)
%Combine by AND
[n,q] = size(pV);
kp = ones(n,1);
for i=1:q
    kp = kp .* pV(:,i);
end

function kp = Com_OR(pV)
%Combine by OR
pV = 1 - pV;
kp = 1 - Com_AND(pV);

