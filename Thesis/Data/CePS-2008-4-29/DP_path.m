function [pth, pflg] = DP_path(idx,pgh,F,S,num,sflg)

%find a path by dynamic programming
if nargin<6
    sflg=0;
end
len = length(idx);
if length(find(pgh==idx(end)))>0
    aaa=1;
end
tp.sim = 0;
tp.pt = -1;
for i=1:num
    for j=1:len
        Dm(i,j) = tp;
    end
end
        
Dm(1,1).sim = F(idx(1));
Dm(1,1).pt = idx(1);

%tic;
for i=2:len
    v = idx(i);
    flg = length(find(pgh==v));%see if v is already in the patially built graph
    ig(i) = flg;
    
    if flg>=1
        st = 1;
    else
        st = 2;
    end
    %for k=2:min(num,i)%for every possible length of path
    for k=st:num%for every possible length of path    
    %for k=2:min(num,i)%for every possible length of path
        if flg>=1
            kp = k;
        else
            kp = k - 1;
        end
        sim = 0;
        pt = -1;
        for j=1:i-1%for point before v
            if S(idx(j),v)>0&Dm(kp,j).sim>sim%ajacent
                sim = Dm(kp,j).sim;
                pt = idx(j);
            end
        end
        if pt == -1 
            
        else
            Dm(k,i).sim = sim + F(idx(i));
            Dm(k,i).pt = pt;
        end
    end
end
%disp(num2str(toc));


%find maxsim/length
mx = 0;
pt = -1;
pl = -1;
if sflg==0
    st = 2;
else
    st = 1;
end
for k = st:num
    if Dm(k,len).sim/k>mx&Dm(k,len).pt~=-1
        mx = Dm(k,len).sim/k;
        pl = k;
    end
end

if mx==0
    pflg = -1;
    pth = -1;
    return;
end
pflg = 1;

%find the path (inverse)
pth(1) = idx(end);

px = pl;
py = len;
flg=1;


for i=1:num
    for j=1:len
        sim(i,j) = Dm(i,j).sim;
        ppt(i,j) = Dm(i,j).pt;
    end
end


i=2;
while flg==1
    if px==1%end of path
        flg=0;
    else
        pth(i) = Dm(px,py).pt;
        if ig(py)==0
            px = px -1;
        end
        py = find(idx==pth(i));
        i = i + 1;
    end
    
end
%find the path
%pth(pl) = idx(end);

%px = pl;
%py = len;
%for i = pl:-1:2
%    pth(i-1) = Dm(px,py).pt;
%    if pth(i-1)<0
%        aaa=1;
%    end
%    if ig(py)==0
%        px = px -1;
%    end
%    py = find(idx==pth(i-1));
%end

for i=1:num
    for j=1:len
        sim(i,j) = Dm(i,j).sim;
        ppt(i,j) = Dm(i,j).pt;
    end
end
