clearvars
global A alpha delta eta gamma
%Production parameters
A=1; alpha=0.5; delta=0.1;
%Agent and utility parameters
AGE=16; RET=10; beta=0.80; eta=2; gamma=0.5;
ub=0.1; retb=0.2;
%Discretized grids
kgrid=[0:0.005:0.6]'; ngrid=[0:0.01:0.7]';
smooth=0.8; eps1=0.02; eps2=0.0001;
p_matrix=[0.9 0.1;0.5 0.5];
L_ind=[1,2^1,2^2,2^3,2^4,2^5,2^6,2^7,2^8,2^9,2^10,2^11,2^12,2^13,2^14,2^15,2^16,2^17];
L_max=2^17;
LR_ind=[1 3 7 15 31 63 127 255 511 1023 2047 4095 8191 16383 32767 65535 131071 262143];
%Initialization
k1max_empl=zeros(length(kgrid),AGE);
k1max_unempl=zeros(length(kgrid),AGE);
k1max_ret=zeros(length(kgrid),AGE);
nmax=zeros(length(kgrid),AGE);
ret_empl=zeros(length(kgrid),AGE);
ret_unempl=zeros(length(kgrid),AGE);
LE1=zeros(2^(AGE-1),AGE); kkE1=zeros(2^(AGE-1),AGE);
LU1=zeros(2^(AGE-1),AGE); kkU1=zeros(2^(AGE-1),AGE);
LR1=zeros(2^(AGE-1),AGE); kkR1=zeros(2^(AGE-1),AGE);
LE=zeros(2^(AGE-1),AGE);   kkE=zeros(2^(AGE-1),AGE);
LE_new=zeros(2^(AGE-1),AGE);
kk1E=zeros(2^(AGE-1),AGE);  nnE =zeros(2^(AGE-1),AGE);
LU =zeros(2^(AGE-1),AGE);  kkU=zeros(2^(AGE-1),AGE);
kk1U=zeros(2^(AGE-1),AGE);
LR =zeros((2^AGE)-1,AGE);  kkR=zeros((2^AGE)-1,AGE);
kk1R=zeros((2^AGE)-1,AGE);
LRa=zeros((2^AGE)-1,AGE);  LRb=zeros((2^AGE)-1,AGE);
kkEa=zeros((2^AGE)-1,AGE); kkUa=zeros((2^AGE)-1,AGE);
kkRa=zeros((2^AGE)-1,AGE); kkRb=zeros((2^AGE)-1,AGE);
V_empl=zeros(length(kgrid),length(ngrid),2);
V_unempl=zeros(length(kgrid),2);
V_ret=zeros(length(kgrid),1);
loop=1;
%Initial Guesses
dif(loop)=1; K(loop)=3.08; N(loop)=2.3;

while dif(loop)>eps1
    %Value function iterations
    [~,w,r]=ProductionFunction(K(loop),N(loop));
    V1_empl=zeros(length(kgrid),1);
    V1_unempl=zeros(length(kgrid),1);
    V1_ret=zeros(length(kgrid),1);
    for age=AGE:-1:1
        for k=1:length(kgrid)
            for ret=1:2
                for n=1:length(ngrid)
                    c_empl= kgrid(k)*(1+r) + ngrid(n)*w - kgrid ;
                    c_empl= (c_empl).*(c_empl>0);
                    V_empl(:,n,ret)=UtilityFunction(c_empl+eps2,ngrid(n))+beta.*((p_matrix(1,1).*V1_empl+p_matrix(1,2).*V1_unempl).*(ret==1)+(V1_ret).*(ret==2));
                end
                c_unempl= kgrid(k)*(1+r) + ub*w - kgrid ;
                c_unempl= (c_unempl).*(c_unempl>0);
                V_unempl(:,ret)=UtilityFunction(c_unempl+eps2,0)+beta.*((p_matrix(2,1).*V1_empl+p_matrix(2,2).*V1_unempl).*(ret==1)+(V1_ret).*(ret==2));
            end
            c_ret= kgrid(k)*(1+r) + retb*w*(age>=RET) - kgrid ;
            c_ret= (c_ret).*(c_ret>0);
            V_ret(:,1)=UtilityFunction(c_ret+eps2,0)+beta.*V1_ret;
            if age==AGE
            V_empl=V_empl(:,:,1);
            V_unempl=V_unempl(:,1);
            end
            [Vmax_empl(k,1),Ind_V_empl(k,age)] = max(V_empl(:));
            [k1max_empl(k,age),nmax_empl(k,age),ret_empl(k,age)]=ind2sub(size(V_empl),Ind_V_empl(k,age));
            [k1max_unempl(k,age),ret_unempl(k,age)]=find(V_unempl==max(max(V_unempl)));
            Vmax_unempl(k,1)=V_unempl(k1max_unempl(k,age),ret_unempl(k,age));
            k1max_ret(k,age)=find(V_ret==max(V_ret));
            Vmax_ret(k,1)=V_ret(k1max_ret(k,age),1);
        end
        V1_empl=Vmax_empl;
        V1_unempl=Vmax_unempl;
        V1_ret=Vmax_ret;
    end
kkE(1,1)=1;
kkU(1,1)=1;
kkR_early(1,1)=1;
kkR(1,1)=1;
LE(1,1)=0;
LU(1,1)=1;
LR_early(1,1)=0;
LR(1,1)=0;
for age=1:AGE-1
    for i=1:1:L_ind(age)
        kk1E(i,age)=k1max_empl(kkE(i,age),age);
        nnE(i,age) =nmax_empl(kkE(i,age),age);
        kk1U(i,age)=k1max_unempl(kkU(i,age),age);
    end
    for i=1:1:LR_ind(age)
        kk1R(i,age)=k1max_ret(kkR(i,age),age);
    end
    
    for i=1:1:L_ind(age)
        LE(i,age+1) = LE(i,age) * p_matrix(1,1) * (ret_empl(kkE(i,age),age)==1);
        LU(i,age+1) = LU(i,age) * p_matrix(2,2) * (ret_unempl(kkE(i,age),age)==1);
        kkE(i,age+1)=kk1E(i,age);
        kkU(i,age+1)=kk1U(i,age);
        
        LE1(i,age+1) = LU(i,age) * p_matrix(2,1) * (ret_unempl(kkE(i,age),age)==1);
        LU1(i,age+1) = LE(i,age) * p_matrix(1,2) * (ret_empl(kkE(i,age),age)==1);
        kkEa(i,age+1)=kk1U(i,age);
        kkUa(i,age+1)=kk1E(i,age);
        
        LRa(i,age+1) = LU(i,age) * (ret_unempl(kkU(i,age),age)==2) ;
        LRb(i,age+1) = LE(i,age) * (ret_empl(kkE(i,age),age)==2) ;
        kkRa(i,age+1)= kk1U(i,age);
        kkRb(i,age+1)= kk1E(i,age);
    end
    for i=1:1:LR_ind(age)
        LR(i,age+1) = LR(i,age);
        kkR(i,age+1)= kk1R(i,age);
    end
    LE(L_ind(age)+1:L_ind(age+1),age+1) = LE1(1:L_ind(age),age+1);
    LU(L_ind(age)+1:L_ind(age+1),age+1) = LU1(1:L_ind(age),age+1);
    LR(LR_ind(age)+1:LR_ind(age+1),age+1) = [LRa(1:L_ind(age),age+1); LRb(1:L_ind(age),age+1)];
    kkE(L_ind(age)+1:L_ind(age+1),age+1) = kkEa(1:L_ind(age),age+1);
    kkU(L_ind(age)+1:L_ind(age+1),age+1) = kkUa(1:L_ind(age),age+1);
    kkR(LR_ind(age)+1:LR_ind(age+1),age+1) = [kkRa(1:L_ind(age),age+1); kkRb(1:L_ind(age),age+1)];
end
    kkE(kkE==0)=1; kkU(kkU==0)=1; kkR(kkR==0)=1;
    nnE(nnE==0)=1;
    K1(loop)=sum(sum(kgrid(kkE).*LE)+sum(kgrid(kkU).*LU)+sum(kgrid(kkR).*LR));
    N1(loop)=sum(sum(ngrid(nnE).*LE));
    dif(loop+1)=abs(K1(loop)-K(loop));
    loop=loop+1;
    K(loop)=K(loop-1)*smooth+K1(loop-1)*(1-smooth);
    N(loop)=N(loop-1)*smooth+N1(loop-1)*(1-smooth);
end

figure; hold on; plot(sum(LE)); plot(sum(LU)); plot(sum(LR)); 
title('Employment, Unemployment and Retirement Rates by Age');
ylabel('Rates'); ylabel('Age Groups');
legend('Age 15', 'Age 45', 'Age 55');

figure; hold on; plot(nmax(:,15)); plot(nmax(:,45)); plot(nmax(:,55));
xlabel('Asset Holdings');ylabel('Labor Supply Choice');
legend('Age 15', 'Age 45', 'Age 55');

figure; hold on; plot(k1max(:,15)); plot(k1max(:,45)); plot(k1max(:,55));
xlabel('Asset Holdings');ylabel('Savings Choice');
legend('Age 15', 'Age 45', 'Age 55');

function[u]=UtilityFunction(c,n)
global eta gamma
u=(((c.^gamma).*(1-n)^(1-gamma)).^(1-eta)-1)./(1-eta);
end

function[Y,w,r]=ProductionFunction(K,N)
global A alpha delta
Y=A*(K^alpha)*(N^(1-alpha));
w=(1-alpha)*A*(K^alpha)*(N^(-alpha));
r=alpha*A*(K^(alpha-1))*(N^(1-alpha))-delta;
end