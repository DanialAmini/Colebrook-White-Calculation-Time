clear;
clc;

tic

name_=["example4";"example5";"bin";"rome"];
demand_multip=[1,1,0.45,1];

name_n=3;
name_n=4;
%name_n=1;

node=tdfread(strcat(name_(name_n), "\junction.txt")); %read tab delimited file
pipe=tdfread(strcat(name_(name_n), "\pipe.txt"));
reservoir=tdfread(strcat(name_(name_n),"\reservoir.txt"));

f_time('file read done');


tic
node_id=node.ID;
elev=node.Elev; %m
demand=node.Demand/1000; %m^3/s

demand=demand.*demand_multip(name_n);


pipe_id=pipe.ID; 
pipe_node=[pipe.Node1 pipe.Node2]; 
len_=pipe.Length; %m
diam=pipe.Diameter./1000; %m
rough=pipe.Roughness./1000; %m


rsrv_=[reservoir.ID reservoir.Head]; % head(m)


[nr,~]=size(rsrv_);
[np,~]=size(pipe_id);
[nj,~]=size(node_id);

N=nj-nr;



V=0.3048*ones(np,1);



q=ones(np,1);
c1=pi/4;
q=c1*diam.*diam.*V; % m^3/s






nu=0.000001;


h=ones(nj,1);
Re=zeros(np,1);
epsD=zeros(np,1);


epsD=rough./diam;

f_time('initializing variables time');

n_iter=10;
%n_iter=1;
%n_iter=1000;
n_iter=30;
for jj=1:n_iter
    
    fprintf('\n iter no %d\n',jj);
    
    tic
    Re=(4/(pi*nu))*abs(q)./diam;
    
    
    Re_turb=max(Re,4000);
    Re_trans=min(max(Re,2000),4000);
    Re_lam=min(Re,2000);
    
    a_turb=(Re>=4000);
    a_trans=(Re<4000 & Re>2000);
    a_lam=(Re<=2000);
    
    f=zeros(np,1);
    fp=zeros(np,1);
    
    f_time('calculate Re');
    
    res=f_swamee(epsD,Re_turb).*a_turb;
    f=f+res(:,1);
    fp=fp+res(:,2);    
    
    f_time('calculate f swamee');
    
    res=f_dunlop(epsD,Re_trans).*a_trans;
    
    f=f+res(:,1);
    fp=fp+res(:,2);
    
    f_time('calculate f transition');
    
    res=f_lam(epsD,Re_lam).*a_lam;
    f=f+res(:,1);
    fp=fp+res(:,2);
    
    f_time('calculate f laminar');

    %f=f_swamee(epsD,Re);
    %fp=fp_swamee(epsD,Re);
    
    r_turb=(8/(pi^2*9.81)).*len_./diam.^5;
    hL=r_turb.*f.*q.*abs(q);

    qq=abs(q);
    g=r_turb.*qq.*(2.*f+Re.*fp);
    
    f_time('calculate h_L & g ');
    
    tic
    

    nnn=nr+nj; %%%
   
    %A=sparse(nj,nj);
    a=1./g;%%%
    i=pipe_node(:,1); %%
    j=pipe_node(:,2);   %% 
    
    tic
    A=sparse(nnn,nnn);    %%%
    A=A+sparse(i,j,-a,nnn,nnn);  %%%
    A=A+sparse(j,i,-a,nnn,nnn);%%%
    A=A+sparse(i,i,a,nnn,nnn);%%%
    A=A+sparse(j,j,a,nnn,nnn);%%%
    
    A(nj+1:nj+nr,:) = []; %%%
    A(:,nj+1:nj+nr) = [];%%%    
    
    f_time('matrix A');

    
    %F=zeros(nj,1);
    %F=-1.*demand;
    F=zeros(nnn,1); %%%
    a=-1.*demand; %%%
    F(1:nj,1)=a(1:nj,1); %%%
    
    
    b=-q+hL./g;%%%%
    
    for k=1:np
        i=pipe_node(k,1);
        j=pipe_node(k,2);
        F(i)=F(i)+b(k);
        F(j)=F(j)-b(k);
        if i>nj
            F(j)=F(j)+rsrv_(i-nj,2)/g(k);
        elseif j>nj
            F(i)=F(i)+rsrv_(j-nj,2)/g(k);
        end
    end
    
    F(nj+1:nj+nr) = [];%%%    
    
    fprintf('build F, time=%f s\n',toc);

%     p=symrcm(A);
%     %p=symamd(A);
%     A=A(p,p);
%     F=F(p);
%     [~,pp] = sort(p);
%     f_time('matrix reordering');    

    h=A\F;
    f_time('solving sparse system done');
    
%     h=h(pp);
%     f_time('matrix reorder back');

    
    dq=sparse(np,1);
    dq=full(dq);
    for k=1:np
        i=pipe_node(k,1);
        j=pipe_node(k,2);
        a=1/g(k);
        if i<=nj && j<=nj %excluding pipes connected to reservoir
            dq(k)=(hL(k)-h(i)+h(j))/g(k);
        elseif i>nj
            dq(k)=(hL(k)-rsrv_(i-nj,2)+h(j))/g(k);
        elseif j>nj
            dq(k)=(hL(k)-h(i)+rsrv_(j-nj,2))/g(k);
        end
    end

    q=q-dq;
    
    fprintf('updatding variables done , time=%f s\n',toc);


end


q=q.*1000;

disp(q(1));
disp(q(2));
disp(q(np));


function res = f_swamee(epsD,Re)
    y1=5.74./Re.^0.9;
    y2=epsD./3.7+y1;
    c1=-2/log(10);
    y3=c1.*log(y2);
    f=1./(y3.*y3);

    y3=c1.*log(y2);
    fp=(2*0.9*c1).*y1.*f./(y2.*y3.*Re);    
    
    res(:,1)=f;
    res(:,2)=fp;
end



function res=f_dunlop(epsD,Re)
    %aa=-1.5634601348517065795;
    aa=-2*0.9*2/log(10);
    %ab=0.00328895476345399058690;
    ab=5.74/4000^0.9;
    c2=-2/log(10);
    
    ac=aa*ab;
    
    y2=epsD./3.7+ab;
    y3=c2.*log(y2);
    
    fa=1./(y3.*y3);
    fb=(2+ac./(y2.*y3)).*fa;
    
    x1=7.*fa-fb;
    x2=0.128-17.*fa+2.5.*fb;
    x3=-0.128+13.*fa-(fb+fb);
    x4=0.032-3.*fa+0.5.*fb;
    
    R=Re./2000;
      
    f=x1+R.*(x2+R.*(x3+R.*x4));
    fp=x2+R.*(2.*x3+3.*R.*x4);
    fp=fp./2000;
    
    res(:,1)=f;
    res(:,2)=fp;
end



function res=f_lam(epsD,Re)
    res(:,1)=64./Re;
    res(:,2)=-64./(Re.*Re);
end


function res=f_time(a)
    res=toc;
    fprintf(strcat(a,', time=%f s\n'),res);
    tic
end


