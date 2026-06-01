function r=BayesianIV_gaussian_fast_flatg2_gprior(Y,X,Z,C,M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes M draws from the posterior of a Bayesian IV model 
% with Gaussian priors on the first-stage regression. The model has the
% form:
%
% X = C*xi + Z*alpha + u = W*pi + u
% Y = X*beta + C*theta + u*delta + eps
%
% Y:        Tx1  vector of dependent variable
% X:        Tx1  matrix of regressors in the 2nd-stage equation
% Z:        Txl matrix of instruments (prior ~ N(0,s2*gam2*I))
% C:        Txk matrix of controls (prior ~ N(0,s2*const*I), const=1e4)
% M:        number of draws from the posterior
% 
% The method uses a flat prior on gamma2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% standardize the data
[z mZ sZ]=zscore(Z,1);
[x mX sX]=zscore(X,1);
[y mY sY]=zscore(Y,1);
if isempty(C); 
    C=ones(length(Y),1); c=C; m1C=0; s1C=1; % adding the constant
else
    [c,mC,sC]=zscore(C,1); c = [ones(length(Y),1) c]; m1C=[0 mC]; s1c=[1 sC];
    C = [ones(length(Y),1) C];              % adding the constant
end 

% compute dimensions and form some useful matrices
[T,k]=size(z);
[T,l]=size(c);
%w=[c z];
P=[x c z y]'*[x c z y];
ww=P(2:end-1,2:end-1);
wx=P(2:end-1,1);
wy=P(2:end-1,end);
wc=P(2:end-1,2:l+1);
xcz2=P(1:end-1,1:end-1);
xczy=P(1:end-1,end);
yy=P(end,end);
cx=P(2:l+1,1);
cz=P(2:l+1,1+l+1:end-1);
cy=P(2:l+1,end);
cc=P(2:l+1,2:l+1);
zx=P(1+l+1:end-1,1);
zz=P(1+l+1:end-1,1+l+1:end-1);
zy=P(1+l+1:end-1,end);
zc=P(1+l+1:end-1,2:l+1);

% matrices to store output 
BETA=NaN(M,1);
DELTA=NaN(M,1);

% starting values for Gibbs
gam2draw=.1;
zzT=(z'*z)/T;
invVprior=zeros(l+k); invVprior(end-k+1:end,end-k+1:end)=zzT;%diag([zeros(l,l);zzT]);
pidraw=(ww+invVprior/gam2draw)\wx;
uu=[1;-pidraw]'*xcz2*[1;-pidraw];

for i=1:M
    
    R=[ [eye(1+l);zeros(k,1+l)] [1;-pidraw] ];
         % invxcu2=inv(R'*xcz2*R+diag([.0001^2 zeros(1,l) .0001^2]));   % tiny regularization for numerical stability
         % bhatols=invxcu2*R'*xczy;
    xcu2=R'*xcz2*R+diag([.0001^2 zeros(1,l) .0001^2]);                  % faster and more stable than indented lines (avoids computing inv(xcu2))
    Lxcu2 = chol(xcu2,'lower');                                         % faster and more stable than indented lines
    bhatols=Lxcu2'\(Lxcu2\R')*xczy;                                     % faster and more stable than indented lines
    ssr=yy-xczy'*R*bhatols;    
    s2epsdraw=1./gamrnd((T-2-2-l)/2,2/ssr);
         %betadeltadraw=(randn(1,l+2)*chol(s2epsdraw*invxcu2))'+bhatols;    
    betadeltadraw = sqrt(s2epsdraw)*(Lxcu2'\randn(l+2,1)) + bhatols;    % faster and more stable than indented lines 
    betadraw=betadeltadraw(1);
    BETA(i)=betadraw;
    deltadraw=betadeltadraw(end);
    DELTA(i)=deltadraw;
        
    % block (ii): s2u
    if k > 2
        s2udraw=1./gamrnd((T+k-2)/2,2/(uu+pidraw'*invVprior*pidraw/gam2draw));    
    else 
        s2udraw=1./gamrnd(T/2,2/uu); 
    end
    
    % version of the algoritm with 4 blocks - less efficiant but faster 
    % block (iii): gam2  
    if k > 2
        gam2draw=1./gamrnd( (k-2)/2 , 2/(pidraw(l+1:end)'*zzT*pidraw(l+1:end)/s2udraw) );
    end 

    % block (iv): pi
    if k > 2
        alpha=(deltadraw^2*s2udraw+s2epsdraw)/s2epsdraw;
        B=wx-(wy-wx*(betadraw+deltadraw)-wc*betadeltadraw(2:end-1))*deltadraw*s2udraw/s2epsdraw;
        A=alpha*ww+invVprior/gam2draw;
            % invA=inv(A); 
            % pidraw=(randn(1,k+l)*chol(s2udraw*invA))'+invA*B;
        LA = chol(A,'lower');                               % faster and more stable than indented lines (avoids computing inv(A))
        pihat = LA'\(LA\B);                                 % faster and more stable than indented lines
        pidraw = sqrt(s2udraw)*(LA'\randn(k+l,1)) + pihat;  % faster and more stable than indented lines 
        uu=[1;-pidraw]'*xcz2*[1;-pidraw];
    else 
        aux=(deltadraw^2*s2udraw+s2epsdraw)/s2epsdraw;
        Bxi=cx-cz*pidraw(l+1:end)-[cy-cx*(betadraw+deltadraw)-cc*betadeltadraw(2:end-1)+cz*pidraw(l+1:end)*deltadraw]*deltadraw*s2udraw/s2epsdraw;
        Axi=aux*cc;
        LAxi = chol(Axi,'lower');                               % faster and more stable than indented lines (avoids computing inv(A))
        xihat = LAxi'\(LAxi\Bxi);                               % faster and more stable than indented lines
        xidraw = sqrt(s2udraw)*(LAxi'\randn(l,1)) + xihat;      % faster and more stable than indented lines 
        pidraw(1:l)=xidraw;

        Bal=zx-zc*pidraw(1:l)-[zy-zx*(betadraw+deltadraw)-zc*betadeltadraw(2:end-1)+zc*pidraw(1:l)*deltadraw]*deltadraw*s2udraw/s2epsdraw;
        Aal=aux*zz;
        LAal = chol(Aal,'lower');                               % faster and more stable than indented lines (avoids computing inv(A))
        alhat = LAal'\(LAal\Bal);                               % faster and more stable than indented lines
        aldraw = sqrt(s2udraw)*(LAal'\randn(k,1)) + alhat;      % faster and more stable than indented lines 
        p=min(1,exp(-(k-2)/2*[log(aldraw'*zz*aldraw) - log(pidraw(l+1:end)'*zz*pidraw(l+1:end))]));
        if rand(1)<p; pidraw(l+1:end)=aldraw; end
        uu=[1;-pidraw]'*xcz2*[1;-pidraw];
    end
end

r.BETA=BETA*sY/sX;    
r.DELTA=DELTA*sY/sX;   
