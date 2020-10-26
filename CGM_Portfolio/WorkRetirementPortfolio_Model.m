clearvars

%!!!!!!!!!!!!!!!!!!!!!!!!
%! Exogenous parameters !
%!!!!!!!!!!!!!!!!!!!!!!!!

%AGE=80+1; AGE0=1; Ret_age_early=42; Ret_age_norm=45; Ret_age_late=50;
eta=10; eps=0.00001; beta=0.96;
%Discretized grids
RS_grid=0:0.01:1; LW_grid=4:404; C_grid=0:0.25:375; N_grid=0:0.2:0.8;
infinity=-1e+10;

%!!!!!!!!!!!!!!!!!!!!!!
%! INVESTMENT RETURNS !
%!!!!!!!!!!!!!!!!!!!!!!
sigma_r=0.157^2;
grid(1,1)= -1.73205080756887;
grid(2,1)=  0.0;
grid(3,1)=  1.73205080756887;
gr=grid*sigma_r^0.5;
gret = 1.06+gr;

RiskyReturns(1,1,:)=[gret(1),gret(2),gret(3)];
Prob=[0.1666666666666 0.6666666666666 0.1666666666666];
RiskFree=1.02;
ret_shock_states=3;

%!!!!!!!!!!!!!!!!
%! LABOR INCOME !
%!!!!!!!!!!!!!!!!

grid(1,1)= -1.73205080756887;
grid(2,1)=  0.0;
grid(3,1)=  1.73205080756887;
sigt_y=0.0738;
eyt=grid*sigt_y^0.5;
f_y=zeros(length(eyt),45);
inc_shock_states=3;
Prob2=[1/6,4/6,1/6];
ret_fac=0.68212;
a=-2.170042+2.700381;
b1=0.16818;
b2=-0.0323371/10;
b3=0.0019704/100;
for ind1=21:65
   avg = exp(a+b1*ind1+b2*ind1^2+b3*ind1^3);
   Income(:,ind1-20) = avg*exp(eyt(:,1));
end
Ret_income= ret_fac.*avg;

%!!!!!!!!!!!!!!!!!!!
%! TERMINAL PERIOD: Age 101
%!!!!!!!!!!!!!!!!!!!

VR1(:,1)=Utility(LW_grid,0,1,eta);
opt_cons_VR(:,81)=LW_grid;
opt_risky_VR(:,81)=0;

%!!!!!!!!!!!!!!!!!!!!!!
%! RETIREMENT PERIODS !
%!!!!!!!!!!!!!!!!!!!!!!

for age=99:-1:65
    %for ind1=2:length(RS_grid)
    %    RS_grid(ind1)=RS_grid(ind1)-rand(1)/1e+4;
    %end
    for LW0=1:length(LW_grid)
            CurrentWealth=LW_grid(LW0);
            C_grid_new=C_grid;
            RS_grid_new=RS_grid;
            Investment=repmat(((LW_grid(LW0)-C_grid_new)'),1,length(RS_grid_new),ret_shock_states); Investment=Investment.*(Investment>=0);
            C_grid_new=C_grid_new.*((LW_grid(LW0)-C_grid_new)>=0);
            Returns=repmat(RS_grid_new,length(C_grid_new'),1,ret_shock_states).*repmat(RiskyReturns,length(C_grid_new),length(RS_grid_new),1)...
                +(1-repmat(RS_grid_new,length(C_grid_new'),1,ret_shock_states)).*repmat(RiskFree,length(C_grid_new),length(RS_grid_new),ret_shock_states);
            FutureWealth=Investment.*Returns + repmat(Ret_income,length(C_grid_new'),length(RS_grid_new),ret_shock_states);
            FutureWealth=max(FutureWealth,min(LW_grid)); FutureWealth=min(FutureWealth,max(LW_grid));
            VR1_1=Prob(1)*interp1(LW_grid',VR1,FutureWealth(:,:,1),'spline');
            VR1_2=Prob(2)*interp1(LW_grid',VR1,FutureWealth(:,:,2),'spline');
            VR1_3=Prob(3)*interp1(LW_grid',VR1,FutureWealth(:,:,3),'spline');
            VR=repmat(Utility(C_grid_new',0,1,eta),1,length(RS_grid_new)) + beta.* ( VR1_1 + VR1_2 + VR1_3 );
            VR = max(VR,infinity);
            [ind_cons,ind_risky]=find(max(VR(:))==VR,1,'first');
            opt_cons_VR(LW0,age)=C_grid_new(ind_cons);
            opt_risky_VR(LW0,age)=RS_grid_new(ind_risky);
            opt_Invest_VR(LW0,age)=LW_grid(LW0)-opt_cons_VR(LW0,age);
            maxVR(LW0,1)=max(VR(:));
    end
    VR1(:,1)=maxVR;
%     age %Uncomment this line to see the progress as the algorithm progresses
end

%!!!!!!!!!!!!!!!!!!!
%! WORKING PERIODS !
%!!!!!!!!!!!!!!!!!!!

for age=64:-1:20
    for LW0=1:length(LW_grid)
            CurrentWealth=LW_grid(LW0);
            C_grid_new=C_grid;
            RS_grid_new=RS_grid;
            Investment=repmat((LW_grid(LW0)-C_grid_new)',1,length(RS_grid_new),ret_shock_states,inc_shock_states);
            Returns=repmat(RS_grid_new,length(C_grid_new'),1,ret_shock_states,inc_shock_states).*repmat(RiskyReturns,length(C_grid_new),length(RS_grid_new),1,inc_shock_states)...
                +(1-repmat(RS_grid_new,length(C_grid_new'),1,ret_shock_states,inc_shock_states)).*repmat(RiskFree,length(C_grid_new),length(RS_grid_new),ret_shock_states,inc_shock_states);
            FutureWealth=Investment.*Returns + repmat(reshape(Income(:,age+1-20),1,1,1,inc_shock_states),length(C_grid_new'),length(RS_grid_new),ret_shock_states,1);
            FutureWealth=max(FutureWealth,min(LW_grid)); FutureWealth=min(FutureWealth,max(LW_grid));
            VR1_11=Prob(1)*Prob2(1)*interp1(LW_grid',VR1,FutureWealth(:,:,1,1),'spline');
            VR1_12=Prob(1)*Prob2(2)*interp1(LW_grid',VR1,FutureWealth(:,:,1,2),'spline');
            VR1_13=Prob(1)*Prob2(3)*interp1(LW_grid',VR1,FutureWealth(:,:,1,3),'spline');
            VR1_21=Prob(2)*Prob2(1)*interp1(LW_grid',VR1,FutureWealth(:,:,2,1),'spline');
            VR1_22=Prob(2)*Prob2(2)*interp1(LW_grid',VR1,FutureWealth(:,:,2,2),'spline');
            VR1_23=Prob(2)*Prob2(3)*interp1(LW_grid',VR1,FutureWealth(:,:,2,3),'spline');
            VR1_31=Prob(3)*Prob2(1)*interp1(LW_grid',VR1,FutureWealth(:,:,3,1),'spline');
            VR1_32=Prob(3)*Prob2(2)*interp1(LW_grid',VR1,FutureWealth(:,:,3,2),'spline');
            VR1_33=Prob(3)*Prob2(3)*interp1(LW_grid',VR1,FutureWealth(:,:,3,3),'spline');
            VR=repmat(Utility(C_grid_new',0,1,eta),1,length(RS_grid_new)) + beta.* ( VR1_11+VR1_12+VR1_13+VR1_21+VR1_22+VR1_23+VR1_31+VR1_32+VR1_33 );
            [opt_cons_VR(LW0,age),opt_risky_VR(LW0,age)]=find(max(max(VR))==VR,1,'first');
            opt_cons_VR(LW0,age)=C_grid_new(opt_cons_VR(LW0,age));
            opt_risky_VR(LW0,age)=RS_grid_new(opt_risky_VR(LW0,age));
            opt_Invest_VR(LW0,age)=LW_grid(LW0)-opt_cons_VR(LW0,age);
            maxVR(LW0,1)=max(VR(:));
    end
    VR1(:,1)=maxVR;
%     age %Uncomment this line to see the progress as the algorithm progresses
end

%!!!!!!!!!!!!!!
%! Figure 2-B !
%!!!!!!!!!!!!!!

figure; hold on;
plot(1:400,opt_risky_VR(1:400,20)); plot(1:400,opt_risky_VR(1:400,30));
plot(1:400,opt_risky_VR(1:400,55)); plot(1:400,opt_risky_VR(1:400,75));
plot(1:400,opt_risky_VR(1:400,99)); xlim([15 350])
title('Risky Share of Portfolio')
legend('Year 20','Year 30','Year 55','Year 75')

%!!!!!!!!!!!!!
%! Functions !
%!!!!!!!!!!!!!

function[u]=Utility(c,n,gamma,eta)
u=((((c).^(gamma)).*(1-n).^(1-gamma)).^(1-eta))./(1-eta); %Last good one
end

%!!!!!!!!!!!!!!!!!!!!!!
%! The End of Program !
%!!!!!!!!!!!!!!!!!!!!!!


%In the case of early retirement, a benefit is reduced 5/9 of one percent 
%for each month before normal retirement age, up to 36 months. If the number
%of months exceeds 36, then the benefit is further reduced 5/12 of one 
%percent per month.
% Ret_avg=0.5; RB_grid=repmat(Ret_avg,1,(Ret_age_late-Ret_age_early+1));
% for i=1:(Ret_age_late-Ret_age_early+1)
%     if i-(Ret_age_norm-Ret_age_early+1)<-3
%         RB_grid(i)=Ret_avg*(1-(5/9)*12*3/100-(Ret_age_norm-Ret_age_early-3-i+1)*(5/12)*12/100);
%     elseif i-(Ret_age_norm-Ret_age_early+1)>=-3 && i-(Ret_age_norm-Ret_age_early+1)<0
%         RB_grid(i)=Ret_avg*(1-(5/9)*12*(Ret_age_norm-Ret_age_early+1-i)/100);
%     elseif i-(Ret_age_norm-Ret_age_early+1)>0
%         RB_grid(i)=Ret_avg*(1+8*(i-(Ret_age_norm-Ret_age_early+1))/100);
%     end
% end

% function [Z,Zprob] = tauchen(N,mu,rho,sigma)
% Z     = zeros(N,1);
% Zprob = zeros(N,N);
% a     = (1-rho)*mu;
% Z(N)  = sqrt(sigma^2 / (1 - rho^2));
% Z(1)  = -Z(N);
% zstep = (Z(N) - Z(1)) / (N - 1);
% 
% for i=2:(N-1)
%     Z(i) = Z(1) + zstep * (i - 1);
% end
% 
% Z = Z + a / (1-rho);
% 
% for j = 1:N
%     for k = 1:N
%         if k == 1
%             Zprob(j,k) = cdf_normal((Z(1) - a - rho * Z(j) + zstep / 2) / sigma);
%         elseif k == N
%             Zprob(j,k) = 1 - cdf_normal((Z(N) - a - rho * Z(j) - zstep / 2) / sigma);
%         else
%             Zprob(j,k) = cdf_normal((Z(k) - a - rho * Z(j) + zstep / 2) / sigma) - ...
%                          cdf_normal((Z(k) - a - rho * Z(j) - zstep / 2) / sigma);
%         end
%     end
% end
% end
% 
% function c = cdf_normal(x)
%     c = 0.5 * erfc(-x/sqrt(2));
% end

% function[Afinal,Bfinal]=gss_bivariate(A1,A4,B1,B4,fun,eps)
% %This function performs a bivariate golden section search method.
% %Inputs are two points for each variable and a bivariate function.
% %For example: [xstar,ystar]=gss_bivariate(-25,40,-100,111,@(x,y)-x^2-y^2,0.00001)
% p=(sqrt(5)-1)/2;
% A2=p*A1+(1-p)*A4;
% B2=p*B1+(1-p)*B4;
% A3=(1-p)*A1+p*A4;
% B3=(1-p)*B1+p*B4;
% f22=fun(A2,B2);
% f23=fun(A2,B3);
% f32=fun(A3,B2);
% f33=fun(A3,B3);
% all_fun=[f22,f23,f32,f33];
% area=10;
% while area>eps
%     if max(all_fun)==f22
%         A1=A1; B1=B1;
%         A4=A3; B4=B3;
%         A3=A2; B3=B2;
%         A2=p*A1+(1-p)*A4;
%         B2=p*B1+(1-p)*B4;
%     elseif max(all_fun)==f33
%         A1=A2; B1=B2;
%         A4=A4; B4=B4;
%         A2=A3; B2=B3;
%         A3=(1-p)*A1+(p)*A4;
%         B3=(1-p)*B1+(p)*B4;
%     elseif max(all_fun)==f23
%         A1=A1; B1=B2;
%         A4=A3; B4=B4;
%         A3=A2; B2=B3;
%         A2=p*A1+(1-p)*A4;
%         B3=(1-p)*B1+(p)*B4;
%     elseif max(all_fun)==f32
%         A1=A2; B1=B1;
%         A4=A4; B4=B3;
%         A2=A3; B3=B2;
%         A3=(1-p)*A1+(p)*A4;
%         B2=p*B1+(1-p)*B4;
%     end
%     f22=fun(A2,B2);
%     f23=fun(A2,B3);
%     f32=fun(A3,B2);
%     f33=fun(A3,B3);
%     all_fun=[f22,f23,f32,f33];
%     area=abs((A4-A1)*(B4-B1));
% end
% Afinal=A2;
% Bfinal=B2;
% end

% function[B]=gss(A,D,fun)
% %This function performs the golden section search method as presented by
% %Algorithm 11.6.1.  Written by Philip Shaw, Fordham University, 2018.
% p=(sqrt(5)-1)/2;
% B=p*A+(1-p)*D;
% C=(1-p)*A+p*D;
% fB=fun(B);
% fC=fun(C);
% dist=2;
% eps=.000001;
% while dist>eps*max([1 abs(B)+abs(C)])
%     if fB>fC
%         D=C;
%         C=B;
%         fC=fB;
%         B=p*C+(1-p)*A;
%         fB=fun(B);
%     else
%         A=B;
%         B=C;
%         fB=fC;  
%         C=p*B+(1-p)*D;
%         fC=fun(C);
%     end
%     dist=abs(D-A);
% end
% end
