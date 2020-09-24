%% Preprocess the spectrum data
%   X_new = Preprocessing(method, X, [...])
%   X: size=nxm (n samples, m frequency points), is the original spectrums data matrix
%   Supported methods and corresponding parameters:
%      Method                                           |  Parameter list
%       'SG'/'Savitzky-Golay'                         	|   halfwidth, power=2, derivative=0
%       'MC'/'mean-centering'                           |   miu=mean(X)
%       'AS'/'auto scaling'/'normal'                    |   miu=mean(X), std=std(X)
%       'SNV'/'standard normal variate'                 |   
%       'RNVx' (robust normal variate, e.g. RNV25, RNV50, RNV75)
%       'AN'/'area normalization'                   	|
%       'VN'/'vector normalization'                     |
%       'MMN'/'maximum minimum normalization'           |           * it's carried out on sample level
%
%       'MSC'/'multiplicative scatter correction'       |   x_ref=mean(X)
%       'OSC'/'orthogonal signal correction'            |   y, ncomp, iter=0, tol=99.9
%       'NAS'/'net analyte signals'                     |   y, ncomp
%
%       'LC'/'linear correction'                        |   index1, index2
%
%       'ALS'/'asymmetric least squares'                |   p=1e-2, lambda=1e7
%       'IALS'/'improved asymmetric least squares'      |   lambda1=1e-3, lambda2=1e3
%
%       'PWT'/'pairwise transfer'                       |   order = 0, x_ref = mean(X)
%       
function [X_output,varargout]=Preprocessing(method,X,varargin)
    [n,m]=size(X);      %   n: sample amount, m: spectrum data dimension
    if ischar(method)
        method=lower(method);
    else
        config=method;
        if strcmp(config.type,'sg')
            [X_output,config]=Preprocessing(method.type,X,config.width,config.power,config.d);
        elseif strcmp(config.type,'mc')
            [X_output,config]=Preprocessing(method.type,X,config.miux);
        elseif strcmp(config.type,'as')
            [X_output,config]=Preprocessing(method.type,X,config.miux,config.stdx);
        elseif strcmp(config.type,'snv')
            [X_output,config]=Preprocessing(method.type,X);
        elseif length(config.type)>3 && strcmp(config.type(1:3),'rnv')
            [X_output,config]=Preprocessing(method.type,X);
        elseif strcmp(config.type,'an')
            [X_output,config]=Preprocessing(method.type,X);
        elseif strcmp(config.type,'vn')
            [X_output,config]=Preprocessing(method.type,X);
        elseif strcmp(config.type,'mmn')
            [X_output,config]=Preprocessing(method.type,X);
        elseif strcmp(config.type,'msc')
            [X_output,config]=Preprocessing(method.type,X,config.x_ref);
        elseif strcmp(config.type,'osc')
            X_centered=X-config.miux;
            X_output=X_centered-X_centered*config.nw*inv(config.np'*config.nw)*config.np'+config.miux; 
        elseif strcmp(config.type,'nas')
            X_output=X*config.F_nap;
        elseif strcmp(config.type,'lc')
            [X_output,config]=Preprocessing(method.type,X,config.index1,config.index2);
        elseif strcmp(config.type,'als')
            [X_output,config]=Preprocessing(method.type,X,config.p,config.lambda);
        elseif strcmp(config.type,'ials')
            [X_output,config]=Preprocessing(method.type,X,config.lambda1,config.lambda2);
        elseif strcmp(config.type,'pwt')
            [X_output,config]=Preprocessing(method.type,X,config.order,config.x_ref);
        else
            error('Unknown preprocessing configuration')
        end
        varargout={config};
        return
    end
    
    if strcmp(method,'savitzky-golay') || strcmp(method,'sg')
        %% Savitzky-Golay smoothing
        width=varargin{1};
        if size(varargin,2)<2
            power=2;
        else
            power=varargin{2};
        end
        if size(varargin,2)<3
            d=0;
        else
            d=varargin{3};
        end
        X_output=SavitzkyGolay(X,width,power,d);
        config.type='sg';
        config.width=width;
        config.power=power;
        config.d=d;
        config.correlated=false;
    end
    %% Normalization, auto scaling, AS
    if strcmp(method,'as')||...
            strcmp(method,'auto scaling')||...
            strcmp(method,'normal')||...
            strcmp(method,'normalization')
        if size(varargin,2)==2
            miux=varargin{1};
            stdx=varargin{2};
        else
            miux=mean(X);
            stdx=std(X);
        end
        X_output=(X-repmat(miux,n,1))./repmat(stdx,n,1);
        config.type='as';
        config.miux=miux;
        config.stdx=stdx;
        config.correlated=true;
    end
    %% Mean centering, MC
    if strcmp(method,'mc')||strcmp('mean centering',method)
        if length(varargin)>0
            miux=varargin{1};
        else
            miux=mean(X);
        end
        X_output=X-miux;
        config.type='mc';
        config.miux=miux;
        config.correlated=false;
    end
    %% Standard Normal Variate transformation, SNV
    if strcmp(method,'snv')||strcmp('standard normal variate',method)
        X_output=(X-repmat(mean(X,2),1,m))./repmat(transpose(std(X')),1,m);
        config.type='snv';
        config.correlated=false;
    end
    %% Robust normal variate, RNVX (e.g. RNV25, RNV50, RNV75)
    if length(method)>3 && strcmp(method(1:3),'rnv')
        percentile_num=str2num(method(4:end));
        X_output=X;
        for i=1:n
            x=X(i,:);
            prct_x=prctile(x,percentile_num);
            std_prct=std(x(x<=prct_x));
            X_output(i,:)=(x-prct_x)/std_prct;
        end
        config.type=method;
        config.correlated=false;
    end
    %% Area normalization, AN
    if strcmp(method,'an')||strcmp(method,'area normalization')
        X_output=X./repmat(sum(X,2),1,m);
        config.type='an';
        config.correlated=false;
    end
    %% Vector Normalization, VN
    if strcmp(method,'vn')||strcmp(method,'vector normalization')
        X=X-repmat(mean(X,2),1,m);    %   zero-centering
        X_output=X./repmat(sqrt(sum(X.^2,2)),1,m);
        config.type='vn';
        config.correlated=false;
    end
    if strcmp(method,'mmn')||strcmp(method,'maximum minimum normalization')
        X_output=X;
        for i=1:n
            x=X(i,:);
            X_output(i,:)=(x-min(x))/(max(x)-min(x));
        end
        config.type='mmn';
        config.correlated=false;
    end
    %% Multiplicative scatter correction, MSC
    if strcmp(method,'msc')
        if size(varargin,2)<1
            x_ref=mean(X);
        else
            x_ref=varargin{1};
        end
        X_output=MSC(X,x_ref);
        config.type='msc';
        config.x_ref=x_ref;
        config.correlated=true;
    end
    %% Orthogonal signal correction, OSC
    if strcmp(method,'osc')||strcmp(method,'orthogonal signal correction')
        y=varargin{1};ncomp=varargin{2};
        if length(varargin)<3
            iter=0;
        else
            iter=varargin{3};
        end
        if length(varargin)<4
            tol=99.9;
        else
            tol=varargin{4};
        end
        addpath('./Dependency/OSC')
        [X_output,nw,np,nt]=osccalc(X-mean(X),y-mean(y),ncomp,iter,tol);
        X_output=X_output+mean(X);
        rmpath('./Dependency/OSC')
        config.type='osc';
        config.ncomp=ncomp;
        config.nw=nw;
        config.np=np;
        config.nt=nt;
        config.miux=mean(X);
        config.stdx=std(X);
        config.correlated=true;
    end
    %% Net analyte signals, NAS
    %   Implemented according to:
    %       Simplification of prediction model for apple sugar content using net analyte preprocessing
    %       张海东,赵杰文,刘木华.利用净分析物预处理法简化苹果糖度预测模型[J].江苏大学学报(自然科学版),2005(04):277-280.
    if strcmp(method,'net analyte signals')||strcmp(method,'nas')||strcmp(method,'nap')
        y=varargin{1};
        ncomp=varargin{2};
        X_neg=[eye(n)-y*((y'*y)^-1)*y']*X;
        [V,D]=eigs(X_neg'*X_neg,ncomp); %   V为特征向量，即U矩阵
        F_nap=eye(m)-V*V';
        X_output=X*F_nap;
        config.type='nas';
        config.F_nap=F_nap;
        config.correlated=true;
    end
    %% Linear
    if strcmp(method,'liniear correction') || strcmp(method,'lc')
        index1=varargin{1};
        index2=varargin{2};
        X_output=LinearCorrection(X,index1,index2);
        config.type='lc';
        config.index1=index1;
        config.index2=index2;
        config.correlated=false;
    end
    %% Asymmetric Least Squares
    if strcmp(method,'asymmetric least squares') || strcmp(method,'als')
        if size(varargin,2)<1
            p=0.01;
        else
            p=varargin{1};
        end
        if size(varargin,2)<2
            lambda=1e7;
        else
            lambda=varargin{2};
        end
        X_output=AsymmetricLeastSquares(X,p,lambda);
        config.type='als';
        config.p=p;
        config.lambda=lambda;
        config.correlated=false;
    end
    %% Improved Asymmetric Least Squares
    if strcmp(method,'improved asymmetric least squares') || strcmp(method,'ials')
        if size(varargin,2)<1
            lambda1=1e-3;
        else
            lambda1=varargin{1};
        end
        if size(varargin,2)<2
            lambda2=1e3;
        else
            lambda2=varargin{2};
        end
        X_output=ImprovedAsymmetricLeastSquares(X,lambda1,lambda2);
        config.type='ials';
        config.lambda1=lambda1;
        config.lambda2=lambda2;
        config.correlated=false;
    end
    %% Pairwise transfer
    if strcmp(method,'pairwise transfer') || strcmp(method,'pwt')
        if length(varargin)>0
            order=varargin{1};
        else
            order=0;
        end
        if length(varargin)>1
            x_ref=varargin{2};
        else
            x_ref=mean(X);
        end
        X_output=pairwise_poly(X,x_ref,order);
        config.type='pwt';
        config.x_ref=x_ref;
        config.order=order;
        config.correlated=true;
    end
    varargout={config};
end
     

%%  SavitzkyGolay平滑
%   对以行组织的等间隔数据y，进行SG平滑，取数据左右各width长度的数据参与运算
%   拟合power阶多项式，并作d阶求导，返回结果
function filterData=SavitzkyGolay(y,width,power,d)

[n,m]=size(y);

xVal=[-width:1:width]';                         %   得到从-w到w的数列
fitMat=[];
for powerID=power:-1:0
    fitMat(:,end+1)=xVal.^powerID;              %   构成如a[-8 -1 0 1 8]+b[4 1 0 1 4]+。。。的形式
end
pinvMat=pinv(fitMat);                           %   求伪逆，如果再右乘数据列向量，就可以得到拟合系数

weightVector=pinvMat(end-d,:)*factorial(d);       %   获得加权系数，是一个行向量

[dataAmount,y_length]=size(y);
filterData=[];

coef=conv(ones(1,m),weightVector);

for i=1:dataAmount
    if d==0
        thisData=conv(y(i,:),weightVector)./coef;   %   这种边缘平滑的技巧只适用于无导数的平滑
    else
        thisData=conv(y(i,:),weightVector);
    end
    filterData(end+1,:)=thisData(width+1:y_length+width);
end
%   为了保险起见，还是不要用这个小聪明了，直接截断吧
filterData=filterData(:,width+1:end-width);
end

%%  ALS非对称最小二乘
%   p：非对称权重，一般取小值，如1e-3~1e-1
%   lambda：光滑性正则，一般较大，如1e2~1e9
%   *   R语言版本https://rdrr.io/cran/baseline/man/baseline.als.html
%   *   默认p=0.01，lambda=7，maxit=20
%   Asymmetric least squares, ALS
%   p: asymmetric weight, a small value, i.e. 1e-3~1e-1
%   lambda: regrularization for smoothness, typically be a large number, i.e. 1e2~1e9
%   *   script implemented by reformatting R script https://rdrr.io/cran/baseline/man/baseline.als.html
function filterData=AsymmetricLeastSquares(y,p,lambda)
    %   首先实现对单行数据的ALS
    epi=1e-5;       %   迭代中止判断条件
    max_step=20;   %   最大迭代次数
    [n,m]=size(y);
    function output=als(y,p,lambda)
        y=y';
        w=ones(1,m);    %   w初始化
        z=mean(y)*ones(m,1);
        D=diff(eye(m),2);
        DD=lambda*D'*D;
        for iter=1:max_step
            W=diag(w);
%             z=pinv(W+DD)*W*y;
            z=(W+DD)\(W*y);
            w_old=w;
            w=p*(y>z)+(1-p)*(y<z);
            if all(abs(w_old-w)<epi)
                break
            end
        end
        output=z';
    end
    filterData=y;
    for i=1:n
        filterData(i,:)=als(y(i,:),p,lambda);
    end
end

%%  IALS改进的非对称最小二乘
%   Improved ALS. (extremely slow!!!)
%   Reference: 姜安,彭江涛,谢启伟,魏纪平,彭思龙.一种改进的非对称最小二乘基线校正算法[J].计算机与应用化学,2012,29(05):537-540.
%   lambda1：一阶导数的平滑参数，一般取小于1e-2
%   lambda2：二阶导数的平滑参数，一般取10~1e4
%   默认lambda1=1e-3，lambda2=1e3，maxit=10
function filterData=ImprovedAsymmetricLeastSquares(y,lambda1,lambda2)
    %   首先实现对单行数据的ALS
    epi=1e-5;       %   迭代中止判断条件
    max_step=20;   %   最大迭代次数
    [n,m]=size(y);
    function output=ials(y,lambda1,lambda2)
        y=y';
        w=ones(1,m);            %   w初始化
        z=mean(y)*ones(m,1);    
        D1=diff(eye(m),1);
        DD1=lambda1*D1'*D1;
        D2=diff(eye(m),2);
        DD2=lambda2*D2'*D2;
        for iter=1:max_step
            W=diag(w);
%             z=pinv(W'*W+DD1+DD2)*(W'*W+DD1)*y;
            z=(W'*W+DD1+DD2)\((W'*W+DD1)*y);
            w_old=w;
            w=1.0*(y<z);
            if all(abs(w_old-w)<epi)
                break
            end
        end
        output=y'-z';
    end
    filterData=y;
    for i=1:n
        filterData(i,:)=ials(y(i,:),lambda1,lambda2);
    end
end

%% 多元散射校正
%   Multiplicative Scatter Correction,MSC
function filterData=MSC(X,x_ref)
    [n,m]=size(X);
    filterData=zeros(n,m);
    for i=1:n
        A=[x_ref',ones(m,1)];
        b=transpose(X(i,:));
        coef=A\b;
        filterData(i,:)=(X(i,:)-coef(2))/coef(1);
    end
end

%%  线性校正(给定两个波数点的索引，将其连线拉至0)
%   Linear correction, pull the line between point index1 and index2 to horizontal axis
function filterData=LinearCorrection(X,index1,index2)
    filterData=zeros(n,m);
    x_range=1:m;
    for i=1:n
        spectra=X(i,:);
        k=(spectra(index2)-spectra(index1))/(index2-index1);
        b=spectra(index1)-index1*k;
        filterData(i,:)=spectra-k*x_range-b;
    end
end

%%  -------------配对拟合变换PWT-----------------------
%   输入：
%       x:      待调整光谱，1行m列
%       x0:     基准光谱，希望x向x0尽量靠拢
%       B:      由基构成的矩阵，m行k列，表示有k个不同的基，例如常数基，一次基，二次基
%   输出：
%       x_hat: 	经过最优匹配后的光谱
%       mincost:匹配产生的误差平方和
function [x_hat,varargout]=pairwise_fit(x,x0,B,varargin)
    x=x';
    x0=x0';
    delta=x-x0;
    H=B'*B;
    f=B'*delta;
    option.Display='off';
    alpha=quadprog(H,f,[],[],[],[],[],[],[],option);
    mincost=0.5*alpha'*H*alpha+f'*alpha+0.5*(delta'*delta);
    x_hat=B*alpha+x;x_hat=x_hat';
    if length(varargin)>0
        figure;hold on
        plot(x0,'k')
        plot(x,'r--')
        plot(x_hat,'r')
        legend('ref','before','after')
    end
    varargout={mincost};
end

%% 以多项式为基进行光谱匹配
%   输入：
%       x:          待调整光谱，1行m列
%       x0:        基准光谱，希望x向x0尽量靠拢
%       order:    匹配阶次，<=0仅调整常数项；1为一次，2为2次，以此类推
%   输出：
%       x_hat:      经过最优匹配后的光谱
%       mincost:  匹配产能的二次误差
function [x_hat,varargout]=pairwise_poly_single(x,x0,order,varargin)
    m=size(x,2);
    seq=(1:m)';
    B=[ones(m,1)];
    for i=1:order
        B=[B,(1e-3)*seq.^i];   %   构建基矩阵
    end
    if isempty(varargin)
        [x_hat,mincost]=pairwise_fit(x,x0,B);
    else
        [x_hat,mincost]=pairwise_fit(x,x0,B,varargin);
    end
    varargout={mincost};
end
        
%% 以多项式为基进行光谱匹配，处理多个样本
%   输入：
%       X:          待调整光谱，n行m列
%       x0:        基准光谱，希望x向x0尽量靠拢
%       order:    匹配阶次，<=0仅调整常数项；1为一次，2为2次，以此类推
%   输出：
%       x_hats:      经过最优匹配后的光谱
%       mincosts:  匹配产能的二次误差
function [x_hats,varargout]=pairwise_poly(X,x0,order)
    [n,m]=size(X);
    x_hats=zeros(n,m);
    mincosts=zeros(n,1);
    for i=1:n
        [x_hat,mincost]=pairwise_poly_single(X(i,:),x0,order);
        x_hats(i,:)=x_hat;
        mincosts(i)=mincost;
    end
    varargout={mincosts};
end