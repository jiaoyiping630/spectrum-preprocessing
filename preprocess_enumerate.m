%   This method tries numerous preprocessing methods
function [fh,varargout]=preprocess_enumerate(X,Y,default_paras)
    warning off
    
    folds=5;                        %   folds for cross-validation
    repetation=50;                  %	repeat time for each configuration
    calib_ratio=0.75;               %   proportion of calibration set
    
    %   Here you define the methods that to be compared
    methods_list={
    'None','SG','SG1D','SG2D','MC','AS','SNV','MMN','RNV25','RNV50','RNV75','AN','VN','MSC','OSC','NAS','ALS',...                     %   Single
    'SG+MC','SG+AS','SG+SNV','SG+MMN','SG+RNV25','SG+RNV50','SG+RNV75','SG+AN','SG+VN','SG+MSC','SG+OSC','SG+NAS','SG+ALS',...                    %   SG+X
    'SG1D+MC','SG1D+AS','SG1D+SNV','SG1D+MMN','SG1D+RNV25','SG1D+RNV50','SG1D+RNV75','SG1D+AN','SG1D+VN','SG1D+MSC','SG1D+OSC','SG1D+NAS','SG1D+ALS',...  %   SG1D+X
    'SG2D+MC','SG2D+AS','SG2D+SNV','SG2D+MMN','SG2D+RNV25','SG2D+RNV50','SG2D+RNV75','SG2D+AN','SG2D+VN','SG2D+MSC','SG2D+OSC','SG2D+NAS','SG2D+ALS',...  %   SG2D+X
    'MC+SG','MC+SG1D','MC+SG2D',...                             %   MC+SG
    'AS+SG','AS+SG1D','AS+SG2D',...
    'SNV+SG','SNV+SG1D','SNV+SG2D',...                          %   SNV+SG
    'MMN+SG','MMN+SG1D','MMN+SG2D',...
    'AN+SG','AN+SG1D','AN+SG2D',...                             %   AN+SGxD (from AN+SG1D)
    'VN+SG','VN+SG1D','VN+SG2D',...                             %   VN+SG
    'MSC+SG','MSC+SG1D','MSC+SG2D',...                          %   MSC+SGxD (from MSC+SG2D)
    'OSC+SG','OSC+SG1D','OSC+SG2D',...                          %   OSC+SGxD (from OSC+SG2D)
    'NAS+SG','NAS+SG1D','NAS+SG2D',...                          %   NAS+SGxD
    'ALS+SG','ALS+SG1D','ALS+SG2D',...                          %   NAS+SGxD
    'OSC+MC','OSC+AS','OSC+SNV','OSC+VN','OSC+AN','OSC+MMN',... %   OSC+normalization
    'VN+MC','VN+SNV','MC+SNV','MC+VN','SNV+VN','SNV+MC',...     %   VN,MC,SNV
    'AN+SG+MC','AN+SG1D+MC','AN+SG2D+MC',...                    %   AN+SGxD+MC (from AN+SG2D+MC)
    'SG+AN+MC','SG1D+AN+MC','SG2D+AN+MC',...                    %   SGxD+AN+MC (from SG2D+AN+MC)
    'MSC+SG+MC','MSC+SG1D+MC','MSC+SG2D+MC',...                 %   MSC+SGxD+MC (from MSC+SG2D+MC)
    'ALS+SNV'
    };

    method_number=length(methods_list);
    
    %   start
    [n,m]=size(X);
    ntrain=floor(n*calib_ratio);
    ntest=n-ntrain;
    
    %%  Single trial
    %   Given sets of X_train,Y_train,X_test,Y_test
    %   this function train PLSR model on training set and returns its RMSEP on test set
    function [rmse_test,err]=get_rmse(X_train,Y_train,X_test,Y_test)
        model=PLSR(X_train,Y_train);
        model=model.cv(folds);
        [rmse_test,err]=model.evaluate(X_test,Y_test);
    end

    %   For each preprocessing configuration, {repetation} times experiments will be done
    %   The RMSEP of each trial will be recorded.
    %   In our paper, following 108 schemes were tested:
    %   ----------------- (none or single method)
    %       None
    %       SG
    %       SG1D
    %       SG2D
    %       MC
    %       AS
    %       SNV
    %       AN
    %       VN
    %       MSC
    %       OSC
    %       NAS
    %       ALS
    %       IALS (removed because too slow)
    %       PWT0
    %       PWT1
    %       PWT2
    %   ----------------- (SG + X)
    %  SG+  MC
    %       AS
    %       SNV
    %       AN
    %       VN
    %       MSC
    %       OSC
    %       NAS
    %       ALS
    %       IALS (removed because too slow)
    %       PWT0
    %       PWT1
    %       PWT2
    %   ----------------- (SG1D + X)
    %  SG1D+MC
    %       AS
    %       SNV
    %       AN
    %       VN
    %       MSC
    %       OSC
    %       NAS
    %       ALS
    %       IALS (removed because too slow)
    %       PWT0
    %       PWT1
    %       PWT2
    %   ----------------- (SG2D + X)
    %  SG2D+MC
    %       AS
    %       SNV
    %       AN
    %       VN
    %       MSC
    %       OSC
    %       NAS
    %       ALS
    %       IALS (removed because too slow)
    %       PWT0
    %       PWT1
    %       PWT2
    %   ----------------- (combination of VN, MC, SNV)
    %       VN+MC
    %       VN+SNV
    %       MC+SNV
    %       MC+VN
    %       SNV+VN
    %       SNV+MC
    %   ----------------- (other relatively rare methods)
    %       ALS+SNV
    %       SNV+SG1D
    %       SG2D+AN+MC
    %       MSC+SG2D+MC
    %       AN+SG1D
    %       AN+SG2D+MC
    %   
    result_table=[];
    err_table={};
    index_table={};
    for method_id=1:method_number
        current_method=methods_list{method_id};
        fprintf(['\nChecking method ',num2str(method_id),'/',num2str(length(methods_list)),' : ',current_method,'  '])
        rmse_list=[];
        err_list={};
        index_list={};
        correlated=true;
        for rep=1:repetation
            rng(rep)
            seq=randperm(ntrain+ntest);
            this_train_idx=seq(1:ntrain);
            this_test_idx=seq(ntrain+1:end);
            if correlated
                %   data shuffling
                this_X_train=X(this_train_idx,:);this_Y_train=Y(this_train_idx,:);
                this_X_test=X(this_test_idx,:);this_Y_test=Y(this_test_idx,:);
                [X_train_processed,X_test_processed,correlated]=...
                preprocessing_compound(current_method,this_X_train,this_Y_train,this_X_test,this_Y_test,default_paras);
                %   cached the processed data to avoid unnecessary computation
                cached_X=[X_train_processed;X_test_processed];  
                cached_X=cached_X(inverse_seq(seq),:);
            else
                %   If the preprocessing scheme is independent on sample level,
                %   the cached result can be reused
                X_train_processed=cached_X(this_train_idx,:);this_Y_train=Y(this_train_idx,:);
                X_test_processed=cached_X(this_test_idx,:);this_Y_test=Y(this_test_idx,:);
            end
            %   evaluation
            [rmse_list(end+1),err_list{end+1}]=get_rmse(X_train_processed,this_Y_train,X_test_processed,this_Y_test);
            index_list{end+1}=this_test_idx;
            fprintf('*')
        end
        result_table(:,end+1)=rmse_list';
        err_table(:,end+1)=err_list';
        index_table(:,end+1)=index_list';
    end
    fprintf('\n')
    fh=figure;
    boxplot(result_table,'Labels',methods_list, 'LabelOrientation', 'inline','boxstyle','filled','PlotStyle','traditional');
    ylabel('RMSEP')
    hold on
    line([0.5,method_number+0.5],[std(Y),std(Y)],'Color','k','LineStyle','--')
    varargout={methods_list,result_table,err_table,index_table};
end

function seq_inv=inverse_seq(seq)
    seq_length=length(seq);
    seq_inv=zeros(1,seq_length);
    for i=1:length(seq)
        seq_inv(seq(i))=i;
    end
end

%%
%   Apply single preprocess method on given data
%   default_paras defines the default parameter for SGxD, OSC and NAS
function [X_train_processed,X_test_processed,correlated]=...
    preprocessing_single(method,X_train,Y_train,X_test,Y_test,default_paras)
    sg_0d_w=default_paras.sg_w;
    sg_0d_p=default_paras.sg_p;
    sg_1d_w=default_paras.sg1d_w;
    sg_1d_p=default_paras.sg1d_p;
    sg_2d_w=default_paras.sg2d_w;
    sg_2d_p=default_paras.sg2d_p;
    osc_comp=default_paras.osc_c;
    nas_comp=default_paras.nas_c;
    if strcmp(method,'None')
        X_train_processed=X_train;
        X_test_processed=X_test;
        correlated=false;
    else
        if strcmp(method,'SG')
            [X_train_processed,config]=Preprocessing('sg',X_train,sg_0d_w,sg_0d_p,0);
        elseif strcmp(method,'SG1D')
            [X_train_processed,config]=Preprocessing('sg',X_train,sg_1d_w,sg_1d_p,1);
        elseif strcmp(method,'SG2D')
            [X_train_processed,config]=Preprocessing('sg',X_train,sg_2d_w,sg_2d_p,2);
        elseif strcmp(method,'MC')
            [X_train_processed,config]=Preprocessing('mc',X_train);
        elseif strcmp(method,'AS')
            [X_train_processed,config]=Preprocessing('as',X_train);
        elseif strcmp(method,'SNV')
            [X_train_processed,config]=Preprocessing('snv',X_train);
        elseif length(method)>3 && strcmp(method(1:3),'RNV')
            [X_train_processed,config]=Preprocessing(method,X_train);
        elseif strcmp(method,'MMN')
            [X_train_processed,config]=Preprocessing('mmn',X_train);
        elseif strcmp(method,'AN')
            [X_train_processed,config]=Preprocessing('an',X_train);
        elseif strcmp(method,'VN')
            [X_train_processed,config]=Preprocessing('vn',X_train);
        elseif strcmp(method,'MSC')
            [X_train_processed,config]=Preprocessing('msc',X_train);
        elseif strcmp(method,'OSC')
            [X_train_processed,config]=Preprocessing('osc',X_train,Y_train,osc_comp);
        elseif strcmp(method,'NAS')
            [X_train_processed,config]=Preprocessing('nas',X_train,Y_train,nas_comp);
        elseif strcmp(method,'ALS')
            [X_train_processed,config]=Preprocessing('als',X_train);
        elseif strcmp(method,'IALS')
            [X_train_processed,config]=Preprocessing('ials',X_train);
        elseif strcmp(method,'PWT0')
            [X_train_processed,config]=Preprocessing('pwt',X_train,0);
        elseif strcmp(method,'PWT1')
            [X_train_processed,config]=Preprocessing('pwt',X_train,1);
        elseif strcmp(method,'PWT2')
            [X_train_processed,config]=Preprocessing('pwt',X_train,2);
        else
            error('Unsupported method.')
        end
        X_test_processed=Preprocessing(config,X_test);
        correlated=config.correlated;
    end
end

%%  
%   Use a sequence of methods to preprocess the spectrum data
%   different components were separated by '+'
%   correlated represents whether this process is independent on sample level,
%   and can reuse if possible for accelaration
function [X_train_processed,X_test_processed,correlated]=...
    preprocessing_compound(methods,X_train,Y_train,X_test,Y_test,default_paras)
    correlated=false;
    method_list=strsplit(methods,'+');
    X_train_processed=X_train;
    X_test_processed=X_test;
    for i=1:length(method_list)
        [X_train_processed,X_test_processed,this_correlated]=...
            preprocessing_single(method_list{i},X_train_processed,Y_train,X_test_processed,Y_test,default_paras);
        correlated=correlated||this_correlated;
    end
end