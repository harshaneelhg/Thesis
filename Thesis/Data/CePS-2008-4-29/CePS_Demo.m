function [Re] = CePS_Demo(qlist,Wau,au_idx,sz,km,para)


%find Center-Piece Subgraph
%2007-8-25
%Wau: adjacent matrix (co-authorship for DBLP dataset)
%au_idx: the name list for the nodes (author names for DBLP)
%qlist: a cell matrix containing author name
%sz: the size (default is 10)
%km: k_softAnd parameter (default is Hard-AND)
%para: parameter to do random walk w/ restart

%----------------------------------------------------------------
% Copyright: 2007
%         Hanghang Tong, Christos Faloutsos
%         All rights reserved.
% Please address questions or comments to: htong (AT) cs (DOT) cmu (DOT)
% edu
%----------------------------------------------------------------
%%paper:
%        Hanghang Tong, Christos Faloutsos
%        Center-piece subgraphs: problem definition and fast solutions. 
%        KDD 2006: 404-413
%----------------------------------------------------------------


%parameter to do rwr
if nargin<6
    para = [10 0.5 2 1 1 0];
end
%soft-K coefficient
if nargin<5
    %use hard AND
    km = length(qlist);
end
if nargin<4
    %size of CePS
    sz = 10;
end

len = length(qlist);
%find query nodes
for i=1:len
    tmp = find_spos(au_idx,qlist{i});
    if tmp==-1
        fprintf('No such Author: %s in the Db\n',qlist{i});
        Re = -1;
        return;
    else
        ql(i) = tmp;
    end
end

%setup Y 
[m,n] = size(Wau);
Y = sparse(m,len);
for i=1:len    
    Y(ql(i),i) = 1;  
    st(i) = ql(i);
end
%Get S
D0 = sum(Wau);
D0 = max(D0,1e-10);
D = sparse(diag(D0.^(-0.5)));
S = D * Wau * D;
for i=1:len
    [F(:,i), realIter] = ppr_i2(S, 0.9, Y(:,i));
    F(:,i) = F(:,i)/sum(F(:,i));
    
end

%get km-more prob (softand)
oF = Kmore_Prob(F,km,1);


%%find index and candidates
[tmp,id] = sort(full(oF),1,'ascend');
cand = find(oF>tmp(max(1,end-min(1000,10*sz))));
%cand = find(oF>=tmp(1));
[id2,I2] = sort(full(F(cand,:)),1,'descend');
[tmp,id] = sort(full(oF(cand,:)),1,'ascend');

%initial partially built graph
pgh = st;


%extract component: find center-piece guy and good connections
num = sz;
vol=0;
v1 = 0;
cnt = 0;
while vol<sz&cnt<=sz    
    pgh = find_a_path2(pgh,oF,S,max(3,min(num,sz-vol)),id2,I2,id,qlist,cand,F,km);
    vol = length(pgh);
    cnt = cnt+1;    
end


sg = pgh2sg(pgh,oF,F,Wau,S,qlist,au_idx);

    

Re.sg = sg;
Re.km = km;
Re.oF = oF;
Re.F = F;
flag = 1;

%plot the subgraph by graphviz
my_draw_dot(Re.sg);

function pgh = find_a_path2(pgh,oF,S,num,id2,I2,id,qlist,cand,F,km)
%find a path
%pick up most important point
mip = -1;
len = length(id);
i = len;
while mip<=0
    if length(find(pgh==cand(id(i))))==0
        mip = cand(id(i));
        pos = id(i);
        %thre = oF(max(1,id(i-20*num)));
    else
        i = i - 1;
    end
end

%for a path from source i to mip
%totally km path
sflg = 0;
tmp = sort(F(pos,:));
tpos = find(F(pos,:)>=tmp(end-km+1));
for j=1:length(tpos)
    i = tpos(j);
    pt = find(I2(:,i)==pos);
    
    idx = cand(I2(1:pt,i));   
    
    
    [pth pflg]= DP_path(idx,pgh,oF,S,num,sflg);
    %add path to the partially built graph
    if pflg>0
        pgh = union(pgh,pth);
        sflg = 1;
        %pgh = setdiff(pgh,mip);
    end
    
end
%pgh = union(pgh,mip);


function sg = pgh2sg(tp,oF,F,Wau,S,qlist,aul)

for i=1:length(tp)
    name{i} = aul(tp(i)).fname;
end
sg.gph = full(Wau(tp,tp)); 
tmp = sg.gph;
tmp(find(tmp>0))=1;
sg.gph2 = tmp;


sg.score = full(oF(tp))';
sg.score = sg.score/max(sg.score);
sg.name = name;
sg.sim1 = full(F' * F);
sF = full(F(tp,:));
ds = full(sF*sF');
for i=1:length(tp)
    for j=1:length(tp)
        cs(i,j) = ds(i,j)/norm(sF(i,:))/norm(sF(j,:));
    end
end
for i=1:length(tp)
    cs(i,i) = 0;
end
sg.sim2 = cs;
%disp(num2str(toc));
cnt = 1;
for i=1:length(tp)-1
    for j=i+1:length(tp)
        if sg.gph(i,j)>0
            %edge weight
            %from i->j
            Ew1 = S(tp(i),tp(j)) * sg.score(j);
            %from j->i
            Ew2 = S(tp(j),tp(i)) * sg.score(i);
            pr =[sg.name{i} ' (' num2str(sg.score(i)) ')- ',num2str(sg.gph(i,j)), '/',num2str(Ew1+Ew2), '  -' sg.name{j} ' (' num2str(sg.score(j)) ')'];
            pair(cnt,1:length(pr)) = pr;
            cnt = cnt + 1;
        end
    end
end
%if isempty(pair)
%else
%    sg.sim3 = pair;
%end

sg.qlist = qlist;
sg.th = 0;
sg.label = name;



function pos = find_spos(aul,str,flg)
if nargin<3
    flg=1;%from begin 
end
len = length(aul);
pos = -1;
if flg==1
    for i=1:len
        if(strcmp(aul(i).fname,str)==1)
            pos = i;
            return;
        end
    end
else%from end to begin
    for i=len:-1:1
        if(strcmp(aul(i).fname,str)==1)
            pos = i;
            return;
        end
    end
end