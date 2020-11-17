
clear;
clc;
close all;

%Optimizations on; this only works in simulink duh

flatten_ = @(A) reshape(A.',1,[])'; 





NN=3;
%NN=1;
NN=3;


NN=1;


NNN=11;


NNN=3164;

NNN=23;
NNN=10;

NNN=round(sqrt(100)+1);
%NNN=round(sqrt(500)+1);
%NNN=round(sqrt(6000)+1);
%NNN=round(sqrt(25000)+1);
%NNN=round(sqrt(40000)+1);
%NNN=round(sqrt(160000)+1);
%NNN=round(sqrt(10e6)+1);


fprintf('matrix size = %d*%d=%d \n',NNN,NNN,NNN^2);


if NN==1
    N_Re=NNN-1;
    N_epsD=NNN-1;
elseif NN==2 %320 mb memory  0.02
    N_Re=floor(sqrt(9978)+1)-1; 
    N_epsD=floor(sqrt(9978)+1)-1;
elseif NN==3 % 28->35 0.07of memory  1120MB
    N_Re=3164-1;
    N_epsD=3164-1;
elseif NN==3.5  %  28-> 44   2560Mb
    N_Re=5000;
    N_epsD=5000;
elseif NN==4
    N_Re=5000-1; %28->
    N_epsD=5000-1;
end

% N_Re=1000;
% N_epsD=5;


geometric_=1; %1=geometric, 1=arithmetic
rand_=0; %0=regular, 1=rand

[Re,epsD]=initialize_(N_Re,N_epsD,geometric_,rand_);

    disp(' ');
    
    
z_err=zeros(size(Re));

fprintf('10x CW-50iter\n');
tic
N_iter=50;
for j=1:10
    z_=zeros(size(Re));
    z_ = z_CW(Re,epsD,N_iter);
end
toc

disp(' ');

tic
disp('calc error of z_');
z_err=1./sqrt(z_)+2.*log10(2.51./(Re.*sqrt(z_))+epsD./3.7);
disp(['max err=',num2str(max(max(abs(z_err))))])
toc
disp(' ');



for N_iter=[3:5,30]
    fprintf('10x CW %d iter \n',N_iter);
    tic
    for i=1:10
        z_temp = z_CW(Re,epsD,N_iter);
    end
    toc
    err_=MaxRelAbsErr(z_,z_temp,' ~~');
    disp(' ');
end

for N_iter=[2,5]
    fprintf('10x CW newton with %d iter\n',N_iter);
    tic
    for i=1:10
        z_temp = z_CW_Newton(Re,epsD,N_iter);
    end
    toc
    err_=MaxRelAbsErr(z_,z_temp,'   ~~');
    disp('  ');
end

for N_iter=[3,5]
    fprintf('10x CW halley with %d iter\n',N_iter);
    tic
    for i=1:10
        z_temp = z_CW_Halley(Re,epsD,N_iter);
    end
    toc
    err_=MaxRelAbsErr(z_,z_temp,'    ');
    disp('  ');
end

z_3iter=z_CW(Re,epsD,3);
z_4iter=z_CW(Re,epsD,4);
z_5iter=z_CW(Re,epsD,5);



disp(' ');
disp(' existing z methods ');
disp(' ');


fprintf('10x serghides\n');
tic
for i=1:10
    z_sergh=func_z_sergh(Re,epsD);
end
toc
MaxRelAbsErr(z_,z_sergh,'   ');

disp(' ');

fprintf('vatan\n');
tic
for i=1:10
    z_vatan=func_z_vatan(Re,epsD);
end
toc
MaxRelAbsErr(z_,z_vatan,' ');


fprintf('goudard\n');
tic
for i=1:10
    z_goudard=func_z_goudard(Re,epsD);
end
toc
MaxRelAbsErr(z_,z_goudard,'  ');

disp(' ');


fprintf('lambert\n');
tic
for i=1:10
    z_lambert=func_z_lambert(Re,epsD);
end
toc
MaxRelAbsErr(z_,z_lambert,'  ');


disp(' ');


fprintf('praks\n');
tic
for i=1:10
    z_praks=func_z_praks(Re,epsD,1);
end
toc
MaxRelAbsErr(z_,z_praks,'  ');


disp(' ');

fprintf('haaland\n');
tic
for i=1:10
    z_SJH=func_z_other(Re,epsD,2);
end
toc
MaxRelAbsErr(z_,z_SJH,'  ');



disp(' ');

fprintf('swamee\n');
tic
for i=1:10
    z_SJH=func_z_other(Re,epsD,1);
end
toc
MaxRelAbsErr(z_,z_SJH,'  ');


disp(' ');







disp(' ');
disp(' z_prime ');
disp(' ');


disp('10x   zp   CW     \n');
tic
for i=1:10
    zp_=func_zp_(Re,z_);
end
toc


disp(' ');

disp('  zp errors Giu vs. exact zp_ \n'); 
tic
for i=1:10
    zp_giu=func_zp_giu(Re,epsD,z_);
end
toc
MaxRelAbsErr(zp_,zp_giu,'   ');

disp(' ');

fprintf('10x zp swamee \n'); 
tic
for i=1:10
    zp_swamee=func_zp_swamee(Re,epsD,z_,2); 
end
toc
MaxRelAbsErr(zp_,zp_swamee,'     ');


%function_plotter_(Re,epsD,z_,zp_)


lambertx=(log(10)/5.02).*Re;


aa=[flatten_(epsD),flatten_(Re),flatten_(z_),flatten_(zp_)];





function y1 = fu(x)
y1=x;
end


function [Re,epsD] = initialize_(N_Re,N_epsD,geometric_,rand_)

    tic
    disp('create Re & epsD vectors');

    if rand_==0
        x=(0:1/N_Re:1);
        y=(0:1/N_epsD:1);
    else
        x=rand(1,N_Re+1);
        y=rand(1,N_epsD+1);
    end

    if geometric_==1
        x=4000.*exp(log(1e8/4000).*x)';  %geometric, power for this situation is efficient as it is
        y=1e-6.*exp(log(0.05/1e-6).*y)'; %geometric, power for this situation is efficient as it is   
    elseif geometric_==0
        x=ones(N_Re+1,1).*4000+(1e8-4000).*x'; %arithmetic
        y=0.05.*y'; %arithmetic
    end
    toc


    [Re,epsD] = meshgrid(x,y);
    
    Re=Re';
    epsD=epsD';
end


function z_ = z_CW(Re,epsD,N_iter)

    z_=-6.880288946433447.*ones(size(Re)); %1/sqrt(0.028)*log(10)/(-2)
    y1=epsD./3.7;
    y2=-2.180158299154324./Re;     %-2/log(10)*2,51
    for i=1:N_iter
        z_=log(y1+y2.*z_);
    end
    z_= 1.325474527619600./(z_.*z_);  %(ln(10)/(-2))^2
  
end


function z_ = z_CW_Newton(Re,epsD,N_iter)

    z_=-6.880288946433447.*ones(size(Re)); %1/sqrt(0.028)*log(10)/(-2)
    y1=epsD./3.7;
    y2=-2.180158299154324./Re;     %-2/log(10)*2,51
    for i=1:N_iter
        z_=z_-(z_-log(y1+y2.*z_))./(1-y2./(y1+y2.*z_));
    end
    z_= 1.325474527619600./(z_.*z_);  %(ln(10)/(-2))^2
  
end


function z_ = z_CW_Halley(Re,epsD,N_iter)

    z_=-6.880288946433447.*ones(size(Re)); %1/sqrt(0.028)*log(10)/(-2)
    y1=epsD./3.7;
    y2=-2.180158299154324./Re;     %-2/log(10)*2,51
    for i=1:N_iter
        f_=z_-log(y1+y2.*z_);
        temp=y2./(y1+y2.*z_);
        fp_=ones(size(z_))-temp;
        fpp_=-temp.*temp;
        
        z_=z_-2.*f_.*fp_./(2.*fp_.*fp_-f_.*fpp_);
    end
    z_= 1.325474527619600./(z_.*z_);  %(ln(10)/(-2))^2
 
end


function z_ = func_z_praks(Re,epsD,N_iter)

    oo=ones(size(Re));
    delta_=74205.5.*oo+1000.*epsD.*Re;
    A=8.*oo+(2/log(10)).*log(16./Re+epsD./3.7);
    B=-74914381.46./(delta_.*delta_);
    C=1391459721232.67./(delta_.*delta_.*delta_);
    %z_=8.*oo-2.*A./(2.*oo-A.*B);
    z_=8.*oo-A-0.5.*A.*A.*B;
    %z_=8.*oo-(6.*A-3.*A.*A.*B)./(6.*oo-6.*A.*B+A.*A.*C);
    z_=(-2/log(10)).*log(2.51./Re.*z_+epsD./3.7);
    z_=(-2/log(10)).*log(2.51./Re.*z_+epsD./3.7);
    z_=(-2/log(10)).*log(2.51./Re.*z_+epsD./3.7);
    
    z_= 1./(z_.*z_);  %(ln(10)/(-2))^2

  
end



function val_ = MaxRelAbsErr(z_,z_hat,text_)

    val_=max(max(abs(z_-z_hat)./z_))*100;
    disp([text_,',err=', num2str(val_),'%']);
    

end


function z_sergh=func_z_sergh(Re,epsD)

    z_sergh=zeros(size(Re)); % 
    y1=epsD./3.7;
    y_y=2.51./Re;
    y2=-0.868588963806504.*log(y1+12./Re);
    y3=-0.868588963806504.*log(y1+y_y.*y2);
    y4=-0.868588963806504.*log(y1+y_y.*y3);
    y5=y3-y2;
    y5=y5.*y5;
    y6=y2-y5./(y4-2.*y3+y2);
    z_sergh=1./(y6.*y6);

end


function z_vatan=func_z_vatan(Re,epsD)

    z_vatan=zeros(size(Re)); % 
    y1=0.124.*Re.*epsD+log(0.4587.*Re);
    y2=exp(log(y1-0.31).*(y1./(y1+0.9633)));
    y3=0.8686.*log(0.4587.*Re./y2);
    z_vatan=1./(y3.*y3);

end

function z_lambert=func_z_lambert(Re,epsD)
    c1=-2/log(10);
    c2=5.02/log(10);
    c3=1/c2;

    z_lambert=c1.*log(c2.*lambertw(c3.*Re)./Re+epsD./3.71);
    z_lambert=1./(z_lambert.*z_lambert);
end


function z_other=func_z_other(Re,epsD,method_)
    if method_==1  %swamee
        y1=log10(epsD./3.7+5.74./exp(log(Re).*0.9));
        z_other=0.25./(y1.*y1);
    elseif method_==2 %haaland
        y1=-1.8.*log10(exp(log(epsD./3.7).*1.11)+6.9./Re);
        z_other=1./(y1.*y1);
    end
end


function z_goudard=func_z_goudard(Re,epsD)

    oo=ones(size(Re));
    a_g=0.868588963806504;
    b_g=epsD./3.7;
    d_g= 0.458682289441045.*Re;
    s_g=b_g.*d_g+log(d_g);
    q_g=exp(log(s_g).*s_g./(s_g+oo));
    g_g=b_g.*d_g+log(d_g./q_g);
    z_g=log(q_g./g_g);
    g_gg=g_g+oo;
    dla_g=z_g.*g_g./ g_gg;
    dcfa_g=dla_g.*(oo+0.5.*z_g./( g_gg.*g_gg+0.33333333333333.*z_g.*(2.*g_g-oo)));
    z_g=a_g.*(log(d_g./q_g)+dcfa_g);
    z_goudard=1./(z_g.*z_g);

end


function zp_swamee=func_zp_swamee(Re,epsD,z_,Number_)

    zp_swamee=-11.895.*(Re.^-1.9).*(z_.*sqrt(z_)).*exp(1.151./sqrt(z_));

end


function zp_=func_zp_(Re,z_)

    zp_=1./(1./((-2.51*4/log(10)).*z_./(Re.*Re).*exp((log(10)/2)./sqrt(z_)))-0.5.*Re./z_);
    
end


function zp_giu=func_zp_giu(Re,epsD,z_)

    %  zp_ Giustolisi')
    zp_=zeros(size(z_));
    oo=ones(size(z_));
    
    zp_=8.1207.*sqrt(z_)./(epsD.*Re.*sqrt(z_)+9.3492.*oo);
    
    zp_=-2.*zp_./(oo+zp_).*z_./Re;

    zp_giu=zp_;
    
    
end


