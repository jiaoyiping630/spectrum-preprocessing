%%  
%   This function find the best parameter for SG-smoothing
%   Input:
%       X: n x m spectrum data, where n is the sample number and m is the wavenumber number
%       Y: n x 1 property data
%   Output:
%       fh: the figure handle
%       rmse_record: cells of np x nd, 
%           where np is the enumeration of polinomial order, nd is the enumeration of derivative
%           each element is a matix of size repetation x nw, where nw is the enumeration of halfwidth

function [fh,varargout]=sg_enumerate(X,Y)
    warning off
    
    folds=5;                %   number of folds in cross-validation
    repetation=50;          %   repeat multiple trials for robustness
    train_ratio=0.75;       %   the ratio of calibration set in each trial

    halfwidth_list=0:1:12;  %   the half window length, 0:1:12 -> width = 3:2:25
    derivative_list=[0,1,2];%   order of derivative
    power_list=[0,1,2,3];   %   order of polynomial, p = 0 && w = 0 is equivalent to raw spectrum
    
    [n,m]=size(X);
    ntrain=floor(n*train_ratio);
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
    %   Given a specific combination of w,p,d,
    %   this function return mean RMSEP among multiple trials
    function [rmse_mean,rmse_std,rmse_list]=get_avg_rmse(width,power,derivative)
        rmse_list=[];
        for rep=1:repetation
            %   data shuffling
            rng(rep)
            seq=randperm(ntrain+ntest);
            this_train_idx=seq(1:ntrain);
            this_test_idx=seq(ntrain+1:end);
            this_X_train=X(this_train_idx,:);this_Y_train=Y(this_train_idx,:);
            this_X_test=X(this_test_idx,:);this_Y_test=Y(this_test_idx,:);
            %   preprocessing
            this_X_train= Preprocessing('sg',this_X_train,width,power,derivative);
            this_X_test=Preprocessing('sg',this_X_test,width,power,derivative);
            %   build, tune, and evaluate the model
            rmse_test=get_rmse(this_X_train,this_Y_train,this_X_test,this_Y_test);
            rmse_list(end+1)=rmse_test;
            fprintf('*')
        end
        rmse_mean=mean(rmse_list);
        rmse_std=std(rmse_list);
    end

    %%
    rmse_record=cell(length(power_list),length(derivative_list));
    fh=figure;
    for p_id=1:length(power_list)
        for d_id=1:length(derivative_list)
            %   d can not exceed p
            if derivative_list(d_id)>power_list(p_id)
                continue
            end
            subplot(length(power_list),length(derivative_list),(p_id-1)*length(derivative_list)+d_id)
            rmse_w_mean=[];
            rmse_w_std=[];
            %   filter out cases that w is too small
            %   because with 2w + 1 data, we can fit a 2*w order polynomial at most
            valid_width_list=halfwidth_list;
            invalid_flag=(2*halfwidth_list+1<=power_list(p_id));
            valid_width_list(invalid_flag)=[];
            %   rmse evaluation
            rmse_current_record=[];
            for w=valid_width_list
                fprintf('\nChecking parameter w = %d, p = %d, d = %d  ',w,power_list(p_id),derivative_list(d_id))
                [rmse_w_mean(end+1),rmse_w_std(end+1),rmse_record_w]=...
                    get_avg_rmse(w,power_list(p_id),derivative_list(d_id));
                rmse_current_record(end+1,:)=rmse_record_w;
            end
            
            rmse_current_record=rmse_current_record';   %   it would be rep x nw
            rmse_record{p_id,d_id}=rmse_current_record;
            
            %% Plot
            plot(2*valid_width_list+1,rmse_w_mean,'b')
            [minval,idx]=min(rmse_w_mean);idx=idx(1);
            hold on
            plot(2*valid_width_list+1,rmse_w_mean+rmse_w_std,'r--')
            plot(2*valid_width_list+1,rmse_w_mean-rmse_w_std,'r--')
            scatter(2*valid_width_list(idx)+1,rmse_w_mean(idx),'bo')
            xlabel('Window Size')
            ylabel('RMSE')
            title(['p = ',num2str(power_list(p_id)),', d = ',num2str(derivative_list(d_id)),', rmse* = ',num2str(minval)])
            drawnow
        end
    end
    varargout={rmse_record};
    fprintf('\n')
end