


disp(' ');disp('x+y')
for i=1:10
oo=ones(size(z_));

% x=rand(size(z_));
% y=rand(size(z_));

x1=rand(size(z_));
x2=rand(size(z_));
x3=rand(size(z_));
x4=rand(size(z_));

z1=0.0001.*x1;

y1=10.*x1;
y2=10.*x2;
y3=10.*x3;
y4=10.*x4;

xv=flatten_(x1);
%xx=sin(1./x);
yy=log(8.*x1).*x2;

%yy=log(x)+exp(3.*x);

c1=rand();
c2=rand();
c3=rand();
c4=rand();
tic 

%clear x;
%1./sqrt(x);

%clear x;
%size(z_);
%rand(size(z_));
%zeros(size(z_));
%ones(size(z_));
%x=y;
%ones(3164,3164);
%clear x;
%x=y;
%lambertw(c2.*Re);
%x.*sqrt(x);
%x.*y;
%x+y;
%x-y;
%x./y;
%sqrt(x);
%x.*x;
%sqrt(x);
%c1./x;
%x.*x.*x;
%1./(x.*x);

% a=(x.*x);
% a.*a.*x;;

%x.*sqrt(x);
%x./sqrt(x);
%exp(x);
%log(x);
%log10(x);
%exp(log(x).*y);
%exp(c1.*log(x));
%exp(log(c1).*x);
%x.^2.;
%x.^y;
%x.^1.1;
%1.1.^x;
%x.^-0.9;
%-1.*x;
%1./x;
%x.^y;
%x.^4;

%a=x.*x;
%a.*a;

%(c1*c2/log(c3)*exp(c4)).*x;
%(c1*c2/log(c3)*exp(c4))./x;

%(x1.*x2)./x3+exp(x4);
%log(x1)+log(x2)+log(x3);
%log(x1+log(x2+log(x3)));
%exp(-x1-exp(-x2-exp(-x3)));
%exp(x1)+exp(x2)+exp(x3);
%log(0.00001*rand(size(z_)));
%log(1000000.*rand(size(z_)));
%exp(0.00001.*y1);
%exp(log(-8.*x1).*x1);
%log(-8.*x1);

%log(x1)+log(x2)+log(x3)+log(x1+log(x2+log(x3)));
%log(1./x1+log(x2+log(8.*x)));
%x1+x2+x3+y1+y2+y3+y4+x4+xx+yy;
%x1.*x2-x3./y1+y2./(y3+y4-x4).*xx+yy;
%x1.*x2-x3./y1+sqrt(y2)./(y3+y4-x4).*xx+sqrt(yy);

%(c1.*x1+c2.*y1)./(c2.*x2+c1.*y3);
%(c1.*x2+c2.*x2.*x1+oo)./(c3.*x2+c4.*x2.*x2+c1.*oo);

%c1.*oo+c1.*x1+c1.*x2.*x2+c1.*x2.*x2.*x2.*x2+c1.*oo+c2.*x3+c3.*x3.*x3+c4.*x3.*x3.*x3.*x3;

%(c1.*x2+c2./x2.*x1+oo)./(c3.*x2+c4./x2./x2+c1.*oo);

% x1=y2./x1; 
% x2=y1.*x4; 
% x3=x3./x3; 
% y1=x1.*x2; 
% x1=x1.*x4;
% x1=y2./x1; 
% x2=y1.*x4; 
% x3=x3./x3; 
% y1=x1.*x2; 
% x1=x1.*x4;


% x1=(c1.*x2+c2./x2.*x1+oo)./(c3.*x2+c4./x2./x2+c1.*sqrt(x2));
% x2=c1.*oo+c1.*x1+c1.*x2.*x2+c1.*x2.*x2.*x2.*x2 +sqrt(x2).*oo+c2.*x3+c3.*x3.*x3+c4.*x3.*x3.*x3.*x3;
% x3=(c1.*x2+c2./x2.*x1+oo)./(c3.*x2+c4./x2./sqrt(x2)+c1.*oo);
% x2=c1.*oo+c1.*x1+c1.*x2.*x2+c1.*x2.*x2.*sqrt(x2).*x2 +c1.*oo+c2.*x3+c3.*x3.*x3+c4.*x3.*x3.*x3.*x3;
% y1=(c1.*x2+c2./x2.*x1+oo)./(c3.*x2+c4./sqrt(x2)./x2+c1.*oo);
% y3=(c1.*x1+c2.*y1)./(c2*x2+c1*sqrt(x2));

% x1=(c1.*x2+c2./x2.*x1+oo)./(c3.*x2+c4./x2./x2+c1.*(x2));
% x1=c1.*oo+c1.*x1+c1.*x2.*x2+c1.*x2.*x2.*x2.*x2 +(x2).*oo+c2.*x3+c3.*x3.*x3+c4.*x3.*x3.*x3.*x3;
% x1=(c1.*x2+c2./x2.*x1+oo)./(c3.*x2+c4./x2./(x2)+c1.*oo);
% x1=c1.*oo+c1.*x1+c1.*x2.*x2+c1.*x2.*x2.*(x2).*x2 +c1.*oo+c2.*x3+c3.*x3.*x3+c4.*x3.*x3.*x3.*x3;
% x1=(c1.*x2+c2./x2.*x1+oo)./(c3.*x2+c4./(x2)./x2+c1.*oo);
% x1=(c1.*x1+c2.*y1)./(c2*x2+c1*(x2));


%x1=((c1.*x2+c2./x2.*x1+oo)./(c3.*x2+c4./(x2)./x2+c1.*oo))./(c1.*oo+c1.*x1+c1.*x2.*x2+c1.*x2.*x2.*(x2).*x2 +c1.*oo+c2.*x3+c3.*x3.*x3+c4.*x3.*x3.*x3.*x3);

%lambertx+lambertx+lambertx;

%floor(lambertx/1e3)+floor(lambertx/1e5)+floor(lambertx/1e4);
%loglamb=log(lambertx);


%how to vectorize if
%aa=(x1>0.5).*x1;

%floor(y1); %same time as one sum or multiply

%x1=((c1.*x2+c2./x2.*x1+oo)./(c3.*x2+c4./(x2)./x2+c1.*oo))./(c1.*oo+c1.*x1+c1.*x2.*x2+c1.*x2.*x2.*(x2).*x2 +c1.*oo+c2.*x3+c3.*x3.*x3+c4.*x3.*x3.*x3.*x3);
%x1=((c1.*floor(x2)+c2./x2.*x1+oo)./(c3.*x2+c4./(x2)./x2+c1.*oo))./(c1.*oo+c1.*x1+c1.*x2.*x2+c1.*x2.*x2.*(x2).*x2 +c1.*oo+c2.*x3+c3.*x3.*x3+c4.*x3.*x3.*x3.*x3);
%adding one floor function kind of ruins it though
%x1=((c1.*log(x2)+c2./x2.*x1+oo)./(c3.*x2+c4./(x2)./x2+c1.*oo))./(c1.*oo+c1.*x1+c1.*x2.*x2+c1.*x2.*x2.*(x2).*x2 +c1.*oo+c2.*x3+c3.*x3.*x3+c4.*x3.*x3.*x3.*x3);
%x1=((c1.*sqrt(x2)+c2./x2.*x1+oo)./(c3.*x2+c4./(x2)./x2+c1.*oo))./(c1.*oo+c1.*x1+c1.*x2.*x2+c1.*x2.*x2.*(x2).*x2 +c1.*oo+c2.*x3+c3.*x3.*x3+c4.*x3.*x3.*x3.*x3);
%x1=((c1.*(x2<0.5)+c2./x2.*x1+oo)./(c3.*x2+c4./(x2)./x2+c1.*oo))./(c1.*oo+c1.*x1+c1.*x2.*x2+c1.*x2.*x2.*(x2).*x2 +c1.*oo+c2.*x3+c3.*x3.*x3+c4.*x3.*x3.*x3.*x3);
%the < kind of ruins it too but not too much?

%x1=((c1.*(x2<0.5)+c2./x2.*x1+oo)./(c3.*x2+c4./(x2>0.4)./(x2<0.9)+c1.*oo))./(c1.*oo+c1.*x1+c1.*x2.*x2+c1.*x2.*x2.*(x2).*x2 +c1.*oo+c2.*x3+c3.*x3.*x3+c4.*x3.*x3.*x3.*x3);
%<  definitely ruins it

%x2=x2>0.4;
%x1=x1<0.9;
%y1=y2<1;
%x1=((c1.*(x2<0.5)+c2./x2.*x1+oo)./(c3.*x2+c4./(x2>0.4)./(x2<0.9)+c1.*oo))./(c1.*oo+c1.*x1+c1.*x2.*x2+c1.*x2.*x2.*(x2).*x2 +c1.*oo+c2.*x3+c3.*x3.*x3+c4.*x3.*x3.*x3.*x3);
%yes it does mess it up

%x1=((c1.*x2(x2<0.5)+c2./x2.*x1+oo)./(c3.*x2+c4./(x2)./x2+c1.*oo))./(c1.*oo+c1.*x1+c1.*x2.*x2+c1.*x2.*x2.*(x2).*x2 +c1.*oo+c2.*x3+c3.*x3.*x3+c4.*x3.*x3.*x3.*x3);



%xv(xv<0.5); %very slow!

%xv.*(xv<0.5); %much faster
%xv.*(xv<0.5 & xv>0.3); %still fast
%floor(2.*xv);  %still quick? 

%x1=((c1.*x2.^2+c2./x2.*x1+oo)./(c3.*x2.^2+c4./(x2)./x2+c1.*oo))./(c1.*oo.^2+c1.*x1+c1.*x2.*x2+c1.*x2.*x2.*(x2).*x2 +c1.*oo+c2.*x3+c3.*x3.*x3+c4.*x3.*x3.*x3.*x3);
%x1=((c1.*(x2<0.5)+c2./x2.*x1+oo)./(c3.*x2+c4./(x2)./x2+(x2<0.5).*oo))./(c1.*oo+c1.*x1+(x2<0.5).*x2.*(x2<0.5)+c1.*x2.*x2.*(x2).*x2 +c1.*(x2<0.5)+c2.*x3+c3.*x3.*x3+(x2<0.5).*(x2<0.5).*x3.*x3.*x3);
%x1=((floor(y1).*(x2<0.5)+floor(y1)./x2.*x1+oo)./(c3.*x2+c4./(x2)./x2+(x2<0.5).*oo))./(c1.*oo+c1.*floor(y1)+(x2<0.5).*x2.*(x2<0.5)+c1.*x2.*x2.*(x2).*floor(y1) +c1.*(x2<0.5)+c2.*x3+floor(y1).*x3.*x3+(x2<0.5).*(x2<0.5).*x3.*x3.*x3);
%hard to believe that still takes about the time of a single log call
%x1=log(abs(((floor(y1).*(x2<0.5)+floor(y1)./x2.*x1+oo)./(c3.*x2+c4./(x2)./x2+(x2<0.5).*oo))./(c1.*oo+c1.*floor(y1)+(x2<0.5).*x2.*(x2<0.5)+c1.*x2.*x2.*(x2).*floor(y1) +c1.*(x2<0.5)+c2.*x3+floor(y1).*x3.*x3+(x2<0.5).*(x2<0.5).*x3.*x3.*x3)));
%abs takes about 32ms, log takes 100ms

%aa=x1.*(x1>0.5); %very fast
%aa=x1(find(x1>0.3)); %very slow
%aa=x1(x1>0.5); % slow


%now let's test in in a long statement %60ms + 
%x1=((c1.*x2(find(x2>3))+c2./x2.*x1+oo)./(c3.*x2+c4./(x2)./x2+c1.*oo))./(c1.*oo+c1.*x1+c1.*x2.*x2+c1.*x2.*x2.*(x2).*x2 +c1.*oo+c2.*x3+c3.*x3.*x3+c4.*x3.*x3.*x3.*x3);
%not working

%aa=all(xv>0.3);

%aa=x1>x2; fast
sign(x1);

%x1>0.4;

%sqrt(y1);
%log(yy);

%exp(xx);
toc
end



% 
% z_=rand(size(z_));
% disp(' ');disp(' z_.*z_')
% tic
%     for i=1:10
%         z_.*z_; 
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('   aa=epsD.*Re;  ')
% tic
%     for i=1:10
%         epsD.*Re; %17 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('  epsD./Re  ')
% tic
%     for i=1:10
%         epsD./Re; %16 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('  epsD+Re  ')
% tic
%     for i=1:10
%         epsD+Re; %16 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp(' z_.*2   ')
% tic
%     for i=1:10
%         z_.*2; %24 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('2.01.*z_    ')
% tic
%     for i=1:10
%         2.01.*z_; %24 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('  2.01./z_  ')
% tic
%     for i=1:10
%         2.01./z_; %24 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('  z_./2.01  ')
% tic
%     for i=1:10
%         z_./2.01; %24 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('  z_./z_  ')
% tic
%     for i=1:10
%         z_./z_; %24 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('  z_.^1.2  ')
% tic
%     for i=1:10
%         z_.^1.2;  %409 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('  z_.^2  ')
% tic
%     for i=1:10
%         z_.^2;  % 13 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('    z_.*z_ ')
% tic
%     for i=1:10
%         z_.*z_;  % 13 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('  z_.^3  ')
% tic
%     for i=1:10
%         z_.^3;  % 390 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('  aa=z_.*z_; aa=aa.*z_;  ')
% tic
%     for i=1:10
%         aa=z_.*z_;
%         aa=aa.*z_;  % 25 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('  aa=z_.^4  ')
% tic
%     for i=1:10
%         aa=z_.^4;  % 397 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('  aa=z_.*z_; aa=aa.*aa;  ')
% tic
%     for i=1:10
%         aa=z_.*z_;  % 21 ms
%         aa=aa.*aa;
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('  aa=z_.^-2  ')
% tic
%     for i=1:10
%         aa=z_.^-2;  % 397 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp(' 1./(z_.*z_)  ')
% tic
%     for i=1:10
%         aa=1./(z_.*z_);  % 397 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('  log(z_)  ')
% tic
%     for i=1:10
%         log(z_);  % 105 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('  log10(z_)  ')
% tic
%     for i=1:10
%         log10(z_);  % 144 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('   sqrt(z_)  ')
% tic
%     for i=1:10
%         sqrt(z_);  % 32 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('z_.^z_;    ')
% tic
%     for i=1:10
%         z_.^z_;  % 381 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('   exp(log(z_).*z_)  ')
% tic
%     for i=1:10
%         exp(log(z_).*z_);  % 138 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('z_.^1.25;    ')
% tic
%     for i=1:10
%         z_.^1.25;  % 381 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('z_.^1.25;    ')
% tic
%     for i=1:10
%         1.25.^z_;  % 381 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('   exp(log(z_).*z_)  ')
% tic
%     for i=1:10
%         exp(log(z_).*1.25);  % 138 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('   exp(log(z_).*z_)  ')
% tic
%     for i=1:10
%         exp(log(1.25).*z_);  % 138 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('   exp(z_)  ')
% tic
%     for i=1:10
%         exp(z_);  % 138 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('   10.^z_  ')
% tic
%     for i=1:10
%         10.^z_;  % 138 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('   2.3.^z_  ')
% tic
%     for i=1:10
%         2.3.^z_;  % 138 ms
%     end
% toc
% 
% z_=rand(size(z_));
% disp(' ');disp('   max(z_)  ')
% tic
%     for i=1:10
%         max(z_);  % 138 ms
%     end
% toc
% 
% z_=rand(size(z_));
% aa=(log(10)/5.02).*z_;
% disp(' ');disp('   lambertW(z_)  ')
% tic
%     for i=1:10
%         lambertw(aa);  % 138 ms
%     end
% toc
% 
% disp(' ');disp('   z*sqrt(z)  ')
% tic
%     for i=1:10
%         sqrt(z_).*z_;  % 44 ms
%     end
% toc
% 
% disp(' ');disp('   z^1.5 exp(log(z)*1.5)  ')
% tic
%     for i=1:10
%         exp(log(z_).*1.5);  % 165 ms
%     end
% toc
% 
% % 
% % z_=vector of 10mill elements
% % 
% % expression		operation	time in ms
% % z_+z_                 sum         	25
% % z_.*z_                mult            24
% % aa=epsD.*Re;          mult            22
% % epsD./Re              mult            23
% % epsD+Re               sum             23
% % z_.*2                 mult            25
% % 2.01.*z_              mult            24
% % 2.01./z_              div             24
% % z_./2.01              div             24
% % z_./z_                div             25
% % z_.^1.2               powfloat        413
% % z_.^2                 powinteger      25
% % z_.*z_                powfloat        243
% % z_.^3                 powinteger      408
% % aa=z_.*z_;aa=aa.*z_;	multpow         28
% % aa=z_.^4              powinteger      396
% % aa=z_.*z_;aa=aa.*aa;	powmult         21
% % aa=z_.^-2             powintegerneg	395
% % 1./(z_.*z_)           powintegerdiv	13
% % log(z_)               Ln              104
% % log10(z_)             log10           138
% % sqrt(z_)              sqrt            32
% % z_.^z_;               powfloat        381
% % exp(log(z_).*z_)      poww/exp&ln     135
% % z_.^1.25;             powfloat        412
% % z_.^1.25;             powfloat        385
% % exp(log(z_).*z_)      poww/exp&ln     159
% % exp(z_)               exp             48
% % 10.^z_                expbase10       382
% % 2.3.^z_               expbase2.3      381
% % abs(z_)               abs				24
% % max(z_)               max             5
% % lambertw(z_)          lambertW        1202
% % z^1.5                 pow half int    165
% % z*sqrt(z)             using sqrt      44
% 
% 
% 
% 
% 
