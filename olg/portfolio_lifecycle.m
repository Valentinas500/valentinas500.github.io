clearvars
global beta eta gamma eps rf Z2 w Zprob22

%Exogenous parameters
AGE=40; RET=20; eta=10; gamma=0.5; eps=0.00001; beta=0.95;
rf=1.01; er=0.02;
w = spline([1,20,30],[0.5,2,1.5],1:AGE)/10;
%Discretized grids
kgrid=0:0.1:5;
rgrid=0:0.1:1;
stock_shock_states=3;
a1opt=zeros(length(kgrid),length(rgrid),stock_shock_states,AGE);
k1opt=zeros(length(kgrid),length(rgrid),stock_shock_states,AGE);
copt=zeros(length(kgrid),length(rgrid),stock_shock_states,AGE);
kspline=zeros(length(kgrid),length(rgrid),stock_shock_states,AGE);
aspline=zeros(length(kgrid),length(rgrid),stock_shock_states,AGE);
cspline=zeros(length(kgrid),length(rgrid),stock_shock_states,AGE);
V=zeros(length(kgrid),length(rgrid));

%Aproximated density function for returns
[Z2,Zprob2] = tauchen(stock_shock_states,rf+er,0,0.1,1);
for i=1:stock_shock_states
    Z22(1,1,i)=Z2(i); Zprob22(1,1,i)=Zprob2(1,i);
end

V1=zeros(length(kgrid),length(rgrid),stock_shock_states);


for age=AGE:-1:1
    for k0=1:length(kgrid)
        for a0=1:length(rgrid)
            for z0=1:length(Z2)
                c=kgrid(k0)*((1-rgrid(a0))*rf+(rgrid(a0))*Z2(z0)) + w(age) - kgrid; c=c.*(c>0);
                V(:,:)=UtilityFunction(repmat(c',1,length(rgrid)),0) + beta .* sum(repmat(Zprob22,length(kgrid),length(rgrid)) .* V1(:,:,:),3);
                if length(find(V==max(V(:))))>1
                    a1opt(k0,a0,z0,age)=1; [k1opt(k0,a0,z0,age)]=find(V(:,1)==max(V(:,1)));
                    Vmax(k0,a0,z0)=V(k1opt(k0,a0,z0,age),a1opt(k0,a0,z0,age));
                else
                    [k1opt(k0,a0,z0,age),a1opt(k0,a0,z0,age)]=find(V==max(V(:)));
                    Vmax(k0,a0,z0)=V(k1opt(k0,a0,z0,age),a1opt(k0,a0,z0,age));
                end
                copt(k0,a0,z0,age)=kgrid(k0)*(rgrid(a0)*rf+(1-rgrid(a0))*Z2(z0)) + w(age) - kgrid(k1opt(k0,a0,z0,age));
                
                %%%%%%%% If want to use interpolation, uncomment the section below
                if age<40
                        if k1opt(k0,a0,z0,age)<=3
                            kk1=5-4; kk2=5;
                        elseif k1opt(k0,a0,z0,age)>=length(kgrid)-2
                         kk1=length(kgrid)-4; kk2=length(kgrid);
                        else
                            kk1=k1opt(k0,a0,z0,age)-2; kk2=k1opt(k0,a0,z0,age)+2;
                        end
                        if a1opt(k0,a0,z0,age)<=3
                            aa1=5-4; aa2=5;
                        elseif a1opt(k0,a0,z0,age)>=length(rgrid)-2
                            aa1=length(rgrid)-4; aa2=length(rgrid);
                        else
                            aa1=a1opt(k0,a0,z0,age)-2; aa2=a1opt(k0,a0,z0,age)+2;
                        end
                        [kspline(k0,a0,z0,age),aspline(k0,a0,z0,age)]=gss_bivariate(kgrid(kk1),kgrid(kk2),rgrid(aa1),rgrid(aa2),@(kx,ax)interpolate(k0,a0,z0,V1,kx,ax, kgrid, rgrid, age),eps);
                        cspline(k0,a0,z0,age)=kgrid(k0)*(rgrid(a0)*rf+(1-rgrid(a0))*Z2(z0)) + w(age) - kspline(k0,a0,z0,age);
                        Vmaxspline = interpolate(k0,a0,z0,V1,kspline(k0,a0,z0,age),aspline(k0,a0,z0,age),kgrid,rgrid,age);
                end
                %%%%%%%% If want to use interpolation, uncomment the section above
            end
        end
    end
    %V1=Vmax;
    %%%%%%%% If want to use interpolation, uncomment the section below
    if age<40
        V1=Vmaxspline; 
    else
        V1=Vmax;
    end
    %%%%%%%% If want to use interpolation, uncomment the section above
end

%Simulate saving and investment decisions
%for agent who always receives mid-shock
cap1(1)=1;sha1(1)=1;shock1=1;
cap2(1)=1;sha2(1)=1;shock2=2;
cap3(1)=1;sha3(1)=1;shock3=3;
for age=1:1:AGE
cap1(age+1)=k1opt(cap1(age),sha1(age),shock1,age);
sha1(age+1)=a1opt(cap1(age),sha1(age),shock1,age);
con1(age)  =copt(cap1(age),sha1(age),shock1,age);
cap2(age+1)=k1opt(cap2(age),sha2(age),shock2,age);
sha2(age+1)=a1opt(cap2(age),sha2(age),shock2,age);
con2(age)  =copt(cap2(age),sha2(age),shock2,age);
cap3(age+1)=k1opt(cap3(age),sha3(age),shock3,age);
sha3(age+1)=a1opt(cap3(age),sha3(age),shock3,age);
con3(age)  =copt(cap3(age),sha3(age),shock3,age);
end

figure;
subplot(3,3,1); plot(kgrid(cap1));
subplot(3,3,2); plot(rgrid(sha1));
subplot(3,3,3); plot((con1));
subplot(3,3,4); plot(kgrid(cap2));
subplot(3,3,5); plot(rgrid(sha2));
subplot(3,3,6); plot((con2));
subplot(3,3,7); plot(kgrid(cap3));
subplot(3,3,8); plot(rgrid(sha3));
subplot(3,3,9); plot((con3));

function[u]=UtilityFunction(c,p)
global eta eps
u=((c+eps).^(1-eta))./(1-eta)-0.00001*(p>0);
end

function [Z,Zprob] = tauchen(N,mu,rho,sigma,m)
Z     = zeros(N,1);
Zprob = zeros(N,N);
a     = (1-rho)*mu;
Z(N)  = m * sqrt(sigma^2 / (1 - rho^2));
Z(1)  = -Z(N);
zstep = (Z(N) - Z(1)) / (N - 1);

for i=2:(N-1)
    Z(i) = Z(1) + zstep * (i - 1);
end 

Z = Z + a / (1-rho);

for j = 1:N
    for k = 1:N
        if k == 1
            Zprob(j,k) = cdf_normal((Z(1) - a - rho * Z(j) + zstep / 2) / sigma);
        elseif k == N
            Zprob(j,k) = 1 - cdf_normal((Z(N) - a - rho * Z(j) - zstep / 2) / sigma);
        else
            Zprob(j,k) = cdf_normal((Z(k) - a - rho * Z(j) + zstep / 2) / sigma) - ...
                         cdf_normal((Z(k) - a - rho * Z(j) - zstep / 2) / sigma);
        end
    end
end
end

function c = cdf_normal(x)
    c = 0.5 * erfc(-x/sqrt(2));
end    

function[fxhat]=interpolate(k0,a0,z0,V1,kx,ax, kgrid, rgrid, age)
global beta rf Z2 w Zprob22
ck = kgrid(k0)*((1-rgrid(a0))*rf+(rgrid(a0))*Z2(z0)) + w(age) - kx; ck=ck.*(ck>0);
uk = UtilityFunction(ck,0);
fx = interpn(kgrid,rgrid,sum(repmat(Zprob22,length(kgrid),length(rgrid)) .* V1(:,:,:),3),kx,ax,'linear');
fxhat = uk+beta*fx;
end

function[Afinal,Bfinal]=gss_bivariate(A1,A4,B1,B4,fun,eps)
%This function performs a bivariate golden section search method.
%Inputs are two points for each variable and a bivariate function.
%For example: [xstar,ystar]=gss_bivariate(-25,40,-100,111,@(x,y)-x^2-y^2,0.00001)
p=(sqrt(5)-1)/2;
A2=p*A1+(1-p)*A4;
B2=p*B1+(1-p)*B4;
A3=(1-p)*A1+p*A4;
B3=(1-p)*B1+p*B4;
f22=fun(A2,B2);
f23=fun(A2,B3);
f32=fun(A3,B2);
f33=fun(A3,B3);
all_fun=[f22,f23,f32,f33];
area=10;
while area>eps
    if max(all_fun)==f22
        A1=A1; B1=B1;
        A4=A3; B4=B3;
        A3=A2; B3=B2;
        A2=p*A1+(1-p)*A4;
        B2=p*B1+(1-p)*B4;
    elseif max(all_fun)==f33
        A1=A2; B1=B2;
        A4=A4; B4=B4;
        A2=A3; B2=B3;
        A3=(1-p)*A1+(p)*A4;
        B3=(1-p)*B1+(p)*B4;
    elseif max(all_fun)==f23
        A1=A1; B1=B2;
        A4=A3; B4=B4;
        A3=A2; B2=B3;
        A2=p*A1+(1-p)*A4;
        B3=(1-p)*B1+(p)*B4;
    elseif max(all_fun)==f32
        A1=A2; B1=B1;
        A4=A4; B4=B3;
        A2=A3; B3=B2;
        A3=(1-p)*A1+(p)*A4;
        B2=p*B1+(1-p)*B4;
    end
    f22=fun(A2,B2);
    f23=fun(A2,B3);
    f32=fun(A3,B2);
    f33=fun(A3,B3);
    all_fun=[f22,f23,f32,f33];
    area=abs((A4-A1)*(B4-B1));
end
Afinal=A2;
Bfinal=B2;
end