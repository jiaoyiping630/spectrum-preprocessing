%%  
%   This function find the best parameter for OSC and NAS
%   Input:
%       X: n x m spectrum data, where n is the sample number and m is the wavenumber number
%       Y: n x 1 property data
%   Output:
%       fh: the figure handle
%       rmse_record: cells of np x nd, 
%           where np is the enumeration of polinomial order, nd is the enumeration of derivative
%           each element is a matix of size repetation x nw, where nw is the enumeration of halfwidth

function fh=osc_nas_enumerate(X,Y)
    warning off
    
    folds=5;                        %   folds for cross-validation
    repetation=50;                  %	repeat time for each configuration
    calib_ratio=0.75;               %   proportion of calibration set

    max_comp=20;                    %   maximum component will be used (larger comp will cost more time, 
                                    %   whereas litter difference could be observed)
    
    [n,m]=size(X);
    ntrain=floor(n*calib_ratio);
    ntest=n-ntrain;
    
    %%  Single trial
    %   Given sets of X_train,Y_train,X_test,Y_test
    %   this function train PLSR model on training set and returns its RMSEP on test set
    function rmse_test=get_rmse(X_train,Y_train,X_test,Y_test)
        model=PLSR(X_train,Y_train);
        model=model.cv(folds);
        rmse_test=model.evaluate(X_test,Y_test);
    end

    %%
    osc_rmse_all=[];
    nas_rmse_all=[];
    
    for comp=0:max_comp
        fprintf(['\nChecking component = ',num2str(comp),' '])
        osc_rmse_list=[];
        nas_rmse_list=[];
        for rep=1:repetation
            fprintf('*')
            %   shuffle
            rng(rep)
            seq=randperm(ntrain+ntest);
            this_train_idx=seq(1:ntrain);
            this_test_idx=seq(ntrain+1:end);
            this_X_train=X(this_train_idx,:);this_Y_train=Y(this_train_idx,:);
            this_X_test=X(this_test_idx,:);this_Y_test=Y(this_test_idx,:);
            %   for osc
            if comp==0
                this_X_train_osc=this_X_train;
                this_X_test_osc=this_X_test;
            else
                [this_X_train_osc,osc_config]=Preprocessing('OSC',this_X_train,this_Y_train,comp);
                this_X_test_osc=Preprocessing(osc_config,this_X_test);
            end
            osc_rmse_list(end+1)=get_rmse(this_X_train_osc,this_Y_train,this_X_test_osc,this_Y_test);
            %   for nas
            if comp==0
                this_X_train_nas=this_X_train;
                this_X_test_nas=this_X_test;
            else
                [this_X_train_nas,nas_config]=Preprocessing('NAS',this_X_train,this_Y_train,comp);
                this_X_test_nas=Preprocessing(nas_config,this_X_test);
            end
            nas_rmse_list(end+1)=get_rmse(this_X_train_nas,this_Y_train,this_X_test_nas,this_Y_test);
        end
        %   
        osc_rmse_all(:,end+1)=osc_rmse_list';
        nas_rmse_all(:,end+1)=nas_rmse_list';
    end
    fprintf('\n')
    
    %%
    fh=figure;
    subplot(1,2,1);hold on
    plot(0:max_comp,mean(osc_rmse_all),'b')
    plot(0:max_comp,mean(osc_rmse_all)+std(osc_rmse_all),'r--')
    plot(0:max_comp,mean(osc_rmse_all)-std(osc_rmse_all),'r--')
    [min_osc,min_osc_idx]=min(mean(osc_rmse_all));
    scatter(min_osc_idx-1,min_osc,'bo')
    title(['OSC, ncomp* = ',num2str(min_osc_idx-1),', rmse* = ',num2str(min_osc)])
    subplot(1,2,2);hold on
    plot(0:max_comp,mean(nas_rmse_all),'b')
    plot(0:max_comp,mean(nas_rmse_all)+std(nas_rmse_all),'r--')
    plot(0:max_comp,mean(nas_rmse_all)-std(nas_rmse_all),'r--')
    [min_nas,min_nas_idx]=min(mean(nas_rmse_all));
    scatter(min_nas_idx-1,min_nas,'bo')
    title(['NAS, ncomp* = ',num2str(min_nas_idx-1),', rmse* = ',num2str(min_nas)])
end