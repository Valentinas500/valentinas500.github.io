
clearvars

global beta gamma1 eta phi UB SSben ...
    A_M A_P alpha Chi vc zeta1 getta SSben_early d_e d_n ...
    delta A_U A_US gamma2 phiphi

%%%%% ENVIRONMENT
T=18; %Total years possible to live
smoothing=0.80;
mort(1:T)=[0.998562276 0.998023418 0.99808485 0.997567474 0.996324643 0.994532501 0.992382903 ...
    0.989433908 0.983140995 0.972649591 0.952530049 0.926026991 0.887291451 0.819213623 ...
    0.677722672 0.515940907 0.357811933 0.243725901 ]; %Mortality/Survivability

%%%%% PARAMETERS EXOGENOUSLY DETERMINED
% Firm production technology and wages paid
alpha               =0.40; %Production parameter, constant returns to scale
zeta1               =0.50; %Fraction of output paid to workers as wages
% Matching Technology
phi                 =0.50; %Matchning parameter, constant returns to scale
% Agents
beta                =0.96; %Discount factor for agents
eta                 =2.00; %Utility parameter
phiphi              =0.000001;

%%%%% Delta is such that capital-to-output is around 2.5.
delta               =0.45; %Depreciation rate equivalent to 2.085% per year 0.97915^5

%%%%% Gamma1 is such that agents on average spend a third of their time of  
%%%%% or market-related activities.
gamma1              =0.36; %Utility parameter

%%%%% SSben is such that the SS benefits are 40% of average wage
SSben       =0.145;  %SS benefit parameter for retirees.
SSben_early =SSben*0.80; % last: 0.75 (too low),0.80 (too low), 0.85 (too high)

% Retirement disutility
d_e=-0.00;
d_n=-0.00;

%%%%% PARAMETERS such that employment to LF makes sense
% Production technology
A_P                 =1.00; %Production technology
% Matching Technology
A_M                 =1.00; %Job matching technology
% Firm
Chi                 =0.045;  %Separation Rate %Fraction of jobs destroyed each period
vc                  =0.074656; %Cost of posting a vacancy (vc = vacancy cost)
getta               =1.00;  %the increasing cost of vacancy posting
% Agents
gamma2              =2.00; %Utility parameter
A_U                 =1.00; %Utility parameter: increases disutility from work and search if >1
A_US                =1.00; %Utility parameter: increases disutility from search only if >1
UB                  =0.05;
eps2                =5;

%%%%% MAKE INITIAL GUESSES
K=1.4406;
N=2.9518;
S=0.9804;
NExit=0.0457;
n_bar=0.4888;

kgrid              =0:0.01:0.28; 
ngrid              =0:0.1:0.60;       
sgrid              =0:0.1:0.62;

maxVE    =zeros(length(kgrid),length(1:T));k1pointE =zeros(length(kgrid),length(1:T));
npointE  =zeros(length(kgrid),length(1:T));spointE  =zeros(length(kgrid),length(1:T));
retpointE=zeros(length(kgrid),length(1:T));cpointE  =zeros(length(kgrid),length(1:T));
IndE     =zeros(length(kgrid),length(1:T));kkEa=zeros((2^T)-1,T); kkUa=zeros((2^T)-1,T);
maxVU    =zeros(length(kgrid),length(1:T));k1pointU =zeros(length(kgrid),length(1:T));
npointU  =zeros(length(kgrid),length(1:T));spointU  =zeros(length(kgrid),length(1:T));
retpointU=zeros(length(kgrid),length(1:T));cpointU  =zeros(length(kgrid),length(1:T));
IndU     =zeros(length(kgrid),length(1:T));
L_ind=[1,2^1,2^2,2^3,2^4,2^5,2^6,2^7,2^8,2^9,2^10,2^11,2^12,2^13,2^14,2^15,2^16,2^17];
L_max=2^17;
LE1=zeros(2^(T-1),T); kkE1=zeros(2^(T-1),T);LU1=zeros(2^(T-1),T); kkU1=zeros(2^(T-1),T);
LE=zeros(2^(T-1),T);   kkE=zeros(2^(T-1),T);LE_new=zeros(2^(T-1),T);
kk1E=zeros(2^(T-1),T); nnE =zeros(2^(T-1),T);ssE =zeros(2^(T-1),T); rrE =zeros(2^(T-1),T);
ccE =zeros(2^(T-1),T); UtilE=zeros(2^(T-1),T);LU =zeros(2^(T-1),T);  kkU=zeros(2^(T-1),T);
kk1U=zeros(2^(T-1),T); nnU= zeros(2^(T-1),T);ssU= zeros(2^(T-1),T); rrU= zeros(2^(T-1),T);
ccU= zeros(2^(T-1),T); UtilU=zeros(2^(T-1),T);

bigloop            =0; diff_K             =1;
diff_N             =1; diff_S             =1;
diff_NExit         =1; diff_max_DistUnemp =1;
diff_tau_c         =1;

while abs(diff_K)>eps2 || abs(diff_N)>eps2 || abs(diff_S)>eps2 || abs(diff_NExit)>eps2
    bigloop=bigloop+1;

%%%%% INITIALIZE MATRICES FOR HH PROBLEM FOR COMPUTATIONAL SPEED
VE=zeros(length(kgrid),length(ngrid),       1     );
VU=zeros(length(kgrid),     1       ,length(sgrid));
VE1=zeros(length(kgrid),1);
VU1=zeros(length(kgrid),1);

%%%%% SOLVE FOR INTEREST RATE AND WAGE RATE
F=A_P*(K^alpha)*(N^(1-alpha)); %Production
w=zeta1*(1-alpha)*A_P*(K^alpha)*(N^(-alpha)); %Wage
r=alpha*A_P*(K^(alpha-1))*(N^(1-alpha))-delta; %Rental rate of capital

%%%%% SOLVE FIRMS PROBLEM FOR NUMBER OF VACANCIES, V
V=FirmProblemSunday_V2( K, N, S, NExit, n_bar);
Profit=F-w*N-r*K-vc*V;

%%%%% SOLVE HHs' PROBLEM TO GET OPTIMAL POLICIES
for age=T:-1:1
    for k0=1:length(kgrid) %k0=total assets rolled over from the previous period; state
                        for n=1:length(ngrid) %n=labor supply; choice
                            for s=1:length(sgrid) %s=job-search intensity; choice
cE= ( kgrid(k0)*r+ kgrid(k0) ...
    + ngrid(n)*w* ...
    - kgrid ) ;
cE= (cE).*(cE>0);
VE(:,n,1)= Utility(cE,ngrid(n),sgrid(1))' ... 
    + (mort(age)) .* beta.* ( (1-Chi).* (VE1(:,1)) ...
    + (Chi).*(VU1(:,1)));

cU= ( kgrid(k0)*r + kgrid(k0) ...
    + w*UB ...
    - kgrid ) ;
cU= (cU).*(cU>0);
VU(:,1,s)= Utility(cU,ngrid(1),sgrid(s))' ...
    + (mort(age)).* beta.* ...
    ((Matching_Function(sgrid(s), 1, S, V, 1, 2) .*  VE1(:,1)) ...
    + ((1 - Matching_Function(sgrid(s), 1, S, V, 1, 2)).*VU1(:,1))) ;   
                            end
                        end
                        
[maxVE(k0,age),IndE(k0,age)] = max(VE(:));
[k1pointE(k0,age),npointE(k0,age),spointE(k0,age)]=ind2sub(size(VE),IndE(k0,age));
                        
[maxVU(k0,age),IndU(k0,age)] = max(VU(:));
[k1pointU(k0,age),npointU(k0,age),spointU(k0,age)]=ind2sub(size(VU),IndU(k0,age));


    end
    VE1=maxVE(:,age);
    VU1=maxVU(:,age);
end

for age=1:T
    for k0=1:length(kgrid) %k0=total assets rolled over from the previous period; state

cpointE(k0,age)= ( kgrid(k0)*r + kgrid(k0) ... %SIMPLIFIED CAPITAL TAX
    + ngrid(npointE(k0,age))*w ...
    - kgrid(k1pointE(k0,age)) ) ;

cpointU(k0,age)= ( kgrid(k0)*r + kgrid(k0) ...
    + w*UB ...
    - kgrid(k1pointU(k0,age)) ) ;
    end

end


%%%%% POPULATION TRANSITION SIMULATION
%%%%% AND
%%%%% AGENTS' LIVES SIMULATION


%%%%% STARTING CONDITIONS FOR NEWLY BORN
kkE(1,1)=1;
kkU(1,1)=1;
LE(1,1)=0;
LU(1,1)=1;

%%%%% POPULATION TRANSITION SIMULATION FOR YEARS 1 THROUGH T-1
for age=1:T-1
    for i=1:1:L_ind(age)
    %= unemployed
            kk1U(i,age)=k1pointU(kkU(i,age),age);
            ssU(i,age) =spointU(kkU(i,age),age);
            nnU(i,age) =npointU(kkU(i,age),age);
            ccU(i,age) =cpointU(kkU(i,age),age);
            UtilU(i,age)=Utility(ccU(i,age),ngrid(nnU(i,age)),sgrid(ssU(i,age)));
    %= employed
            kk1E(i,age)=k1pointE(kkE(i,age),age);
            ssE(i,age) =spointE(kkE(i,age),age);
            nnE(i,age) =npointE(kkE(i,age),age);
            ccE(i,age) =cpointE(kkE(i,age),age);
            UtilE(i,age)=Utility(ccE(i,age),ngrid(nnE(i,age)),sgrid(ssE(i,age)));
    end

    for i=1:1:L_ind(age)
        LE(i,age+1) = mort(age) * LE(i,age) * (1-Chi) ;
        LU(i,age+1) = mort(age) * LU(i,age) * (1-Matching_Function(sgrid(ssU(i,age)), 1, S, V, 1, 2));
        kkE(i,age+1)=kk1E(i,age);
        kkU(i,age+1)=kk1U(i,age);
        LE1(i,age+1) = mort(age) * LU(i,age) * (Matching_Function(sgrid(ssU(i,age)), 1, S, V, 1, 2)) ;
        LU1(i,age+1) = mort(age) * LE(i,age) * (Chi ) ; 
        kkEa(i,age+1)=kk1U(i,age);
        kkUa(i,age+1)=kk1E(i,age);
    end

    LE(L_ind(age)+1:L_ind(age+1),age+1) = LE1(1:L_ind(age),age+1);
    LE_new(L_ind(age)+1:L_ind(age+1),age+1) = (LE1(1:L_ind(age),age+1)>0);
    LU(L_ind(age)+1:L_ind(age+1),age+1) = LU1(1:L_ind(age),age+1);
    kkE(L_ind(age)+1:L_ind(age+1),age+1) = kkEa(1:L_ind(age),age+1);
    kkU(L_ind(age)+1:L_ind(age+1),age+1) = kkUa(1:L_ind(age),age+1);
end

for age=T
    for i=1:1:L_ind(age)
    %= unemployed
            kk1U(i,age)=k1pointU(kkU(i,age),age);
            ssU(i,age) =spointU(kkU(i,age),age);
            nnU(i,age) =npointU(kkU(i,age),age);
            ccU(i,age) =cpointU(kkU(i,age),age);
            UtilU(i,age)=Utility(ccU(i,age),ngrid(nnU(i,age)),sgrid(ssU(i,age)));
    %= employed
            kk1E(i,age)=k1pointE(kkE(i,age),age);
            ssE(i,age) =spointE(kkE(i,age),age);
            nnE(i,age) =npointE(kkE(i,age),age);
            ccE(i,age) =cpointE(kkE(i,age),age);
            UtilE(i,age)=Utility(ccE(i,age),ngrid(nnE(i,age)),sgrid(ssE(i,age)));
    end
end

kkE(kkE==0)=1; kkU(kkU==0)=1;
nnE(nnE==0)=1; ssU(ssU==0)=1;

%%%%%% COMPUTE NEW AGGREGATE VALUES
%AGGREGATE LABOR
N_pop=sum(LE.*ngrid(nnE));
N1=sum(sum(LE.*ngrid(nnE)));
Empl_age=sum(LE);
Unempl_age=sum(LU);
n_bar1=sum(sum(LE_new.*LE.*ngrid(nnE)))./sum(sum(LE_new.*LE));
%AGGREGATE CAPITAL
K_pop=sum(LE.*kgrid(kkE))+sum(LU.*kgrid(kkU));
K1=sum(sum(LE.*kgrid(kkE)))+sum(sum(LU.*kgrid(kkU)));
%AGGREGATE JOB SEARCH EFFORT
S_pop=sum(LU.*sgrid(ssU));
S1=sum(sum(LU.*sgrid(ssU)));
%LABOR EXIT DUE TO DEATH
NExit1=sum(sum(LE.*(rrE==2).*(ngrid(nnE))))/N;
%NUMBER OF UNEMPLOYED
Omega_pop=(sum(LU.*(ssU>1)));
Omega=sum(sum(LU.*(ssU>1)));
%AGGREGATE CONSUMPTION
CE_pop=sum(LE.*ccE);
CU_pop=sum(LU.*ccU);
C_pop=CE_pop+CU_pop;
C=sum(C_pop);
%POPULATION
P_pop=sum(LE)+sum(LU);
P=sum(P_pop);
K_pop_mean=K_pop./P_pop;
%WELFARE
UT_LE= sum(sum(LE.*UtilE));
UT_LU= sum(sum(LU.*UtilU));
UT= sum(sum(LE.*UtilE))+sum(sum(LU.*UtilU));
UT_LE_pop=(sum(LE.*UtilE));
UT_LU_pop=(sum(LU.*UtilU));
UT_pop=(sum(LE.*UtilE))+   (sum(LU.*UtilU));

%%%%%% COMPUTE DIFFERENCES BETWEEN GUESSES AND COMPUTED AGGREGATE VALUES
diff_K             = (K1/K-1)*100;
diff_N             = (N1/N-1)*100;
diff_S             = (S1/S-1)*100;
diff_n_bar         = (n_bar1/n_bar-1)*100;
diff_NExit         = (NExit1/NExit-1)*100;

%%%%%% STORE AGGREGATE VALUES
store_diff_K(bigloop)             = diff_K;
store_diff_N(bigloop)             = diff_N;
store_diff_S(bigloop)             = diff_S;
store_diff_n_bar(bigloop)         = diff_n_bar;
store_diff_NExit(bigloop)         = diff_NExit;

%%%%%% UPDATE THE INITIAL GUESSES BEFORE REITERATION
K=K*(smoothing)+K1*(1-smoothing);
N=N*(smoothing)+N1*(1-smoothing);
S=S*(smoothing)+S1*(1-smoothing);
n_bar=n_bar*(smoothing)+n_bar1*(1-smoothing);
NExit=NExit*(smoothing)+NExit1*(1-smoothing);

store_K(bigloop)    = K;     store_K1(bigloop)    = K1; 
store_N(bigloop)    = N;     store_N1(bigloop)    = N1;
store_S(bigloop)    = S;     store_S1(bigloop)    = S1;
store_n_bar(bigloop)= n_bar; store_n_bar1(bigloop)= n_bar1;
store_NExit(bigloop)= NExit; store_NExit1(bigloop)= NExit1;
store_V(bigloop)    = V;     store_UT(bigloop)=UT/100;

[store_K;store_N;store_S;store_NExit;store_V]

end

function [U]=Utility(c,n,s)
global gamma1 eta A_U A_US gamma2 phiphi
U=((((c+phiphi).^(gamma1)).*(1-A_U.*(n+A_US.*s.^gamma2)).^(1-gamma1)).^(1-eta)-1)./(1-eta); %Last good one
end

function[RESULT]=Matching_Function(job_search, ~, total_search_effort, V, ~, type)
% Matching function describes the relationship between the unemployed,
% vacancies and job search effort of a particular group of unemployed.
global phi A_M
if type==1
    M= A_M * V^phi * total_search_effort^(1-phi) ;
    q=M/V;
    if q>1
        q=1;%('Prob warning: q>1')
    elseif q<0
        q=0;%('Prob warning: q<1')
    end
    RESULT=q;
elseif type==2
    rel_search_intensity= (job_search) / (total_search_effort);
    M1= rel_search_intensity * A_M * V^phi * total_search_effort^(1-phi) ;
    p=M1/1;
    if p>1
        p=1;%('Prob warning: p>1')
    elseif p<0
        p=0;%('Prob warning: p<1')
    elseif isnan(p)
        p=0;
    end
    RESULT=p;
end
end

function [ V_final ] = FirmProblemSunday_V2( K, N, S, NExit, n_bar)
global phi A_M A_P alpha Chi vc zeta1 delta getta
w=zeta1*(1-alpha)*A_P*(K^alpha)*(N^(-alpha)); %Wage
r=alpha*A_P*(K^(alpha-1))*(N^(1-alpha))-delta; %Rental rate of capital
Vgrid=0.0001:0.0001:15;
FOC=zeros(1,length(Vgrid));
for V=1:length(Vgrid)
    q(V)= A_M * (Vgrid(V)^(phi-1)) * (S^(1-phi));
    N1(V)= (1-Chi-NExit)*N + q(V)*Vgrid(V)*n_bar;
    redpart(V)= (1-alpha) * (A_P*(K^alpha)*(N1(V)^(-alpha))) * q(V) * n_bar ;
    bluepart(V)= q(V)*w*n_bar;
    redblue(V)= (redpart(V)-bluepart(V));
    FOC(V)= ((1+r)^(-1)) * (redpart(V)-bluepart(V)) - getta*vc*V^(getta-1);
end
V_final= Vgrid(abs(FOC)==min(abs(FOC)));
end