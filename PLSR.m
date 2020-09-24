%%  Partial least square regression
%   Parameters:
%       component:  number of component used for regression
classdef PLSR
    properties
        n
        m
        X_train
        Y_train
        component
        coeff
        beta
        
        max_component
        X_mean
        Y_mean
        X0
        Y0
        
        Xloadings
        Yloadings
        Xscores
        Yscores
    end
    methods
        %% Initialization
        function obj=PLSR(X_train,Y_train)
            obj.X_train=X_train;
            obj.Y_train=Y_train;
            [obj.n,obj.m]=size(X_train);
            obj.X_mean=mean(X_train);
            obj.Y_mean=mean(Y_train);
            obj.X0=X_train-repmat(obj.X_mean,obj.n,1);
            obj.Y0=Y_train-repmat(obj.Y_mean,obj.n,1);
            [obj.Xloadings,obj.Yloadings,obj.Xscores,obj.Yscores,...
                BETA,PCTVAR,MSE,stats] = plsregress(obj.X0,obj.Y0);
            obj.coeff=stats.W;
            obj.max_component=size(obj.Xscores,2);  %   这样可直接获得最大主成分数
        end
        %%  Fitting
        function obj=fit(obj,varargin)
            component=varargin{1};
            obj.component=component;
            %   ordinary method, use plsregress each time
            %   but it will be very slow, and has much redundancy 
%             [XL,yl,XS,YS,beta,PCTVAR] = plsregress(obj.X_train,obj.Y_train,component);
%             obj.beta=beta;
            %   calculate a full PLSR first, and the cache information can be reused
            obj.beta=obj.X0'*...
                obj.Yscores(:,1:component)*...
                ((transpose(obj.Xscores(:,1:component))*...
                (obj.X0*obj.X0')*obj.Yscores(:,1:component))\transpose(obj.Xscores(:,1:component))*obj.Y0);
            obj.beta=[obj.Y_mean-obj.X_mean*obj.beta;obj.beta];
        end
        %%  Prediiction
        function prediction=predict(obj,X_test)
            [n,m]=size(X_test);
            prediction=[ones(n,1),X_test]*obj.beta;
        end
        %%  Evaluation
        function [rmse,varargout]=evaluate(obj,X_valid,Y_valid)
            [n,m]=size(X_valid);
            Y_pred=obj.predict(X_valid);
            err=Y_pred-Y_valid;
            rmse=sqrt(sum(sum(err.^2))/n);
            varargout={err,Y_pred};
        end
        %%	Tuning
        function [obj,varargout]=tune(obj,X_valid,Y_valid,varargin)
            
            v1_min_val=1;
            v1_max_val=obj.max_component;
            v1_interval=1;
            v1_name='component';
            v1_range=v1_min_val:v1_interval:v1_max_val;

            result=zeros(1,length(v1_range));
            for i =1:length(v1_range)
                obj=obj.fit(v1_range(i));
                result(i)=obj.evaluate(X_valid,Y_valid);
            end
            [~,k]=min(result);
            bestpara=v1_range(k);
            obj=obj.fit(bestpara);
            varargout={bestpara,result,{v1_name,v1_range}};
        end
        
        %%  Cross-validation
        function [obj,varargout]=cv(obj,folds,varargin)
            
            %   use fixed seed for reproducibility
            store_rng=rng;
            rng(0);
            partition_index=crossvalind('Kfold', obj.n,folds);
            rng(store_rng)
            
            for i=1:folds
                train_idx=(partition_index~=i);
                valid_idx=(partition_index==i);
                model=PLSR(obj.X_train(train_idx,:),obj.Y_train(train_idx,:));
                [~,~,fold_result,parasetting]=model.tune(obj.X_train(valid_idx,:),obj.Y_train(valid_idx,:));
                if i==1
                    result=fold_result;
                else
                    result=utils.add_mismatch(result,fold_result);
                end
            end
            result=result/folds;
            
            if size(varargin,2)>0
                showflag=varargin{1};
            else
                showflag=0;
            end
            
            if length(parasetting)==2
                if length(result)<length(parasetting{2})
                    parasetting{2}=parasetting{2}(1:length(result));
                end
            elseif length(parasetting)==4
                if size(result,1)<parasetting{2}
                    parasetting{2}=parasetting{2}(1:size(result,1));
                end
                if size(result,2)<parasetting{4}
                    parasetting{4}=parasetting{4}(1:size(result,2));
                end
            else
                error('oops! something is wrong ...')
            end
            
            %%  define output
            v1_name=parasetting{1};
            v1_range=parasetting{2};
            [~,folds]=min(result);
            bestpara=v1_range(folds);
            obj=obj.fit(bestpara);
            
            %   build output structure
            out_struct=struct;
            out_struct.variable_names={v1_name};
            out_struct.variable_ranges={v1_range};
            out_struct.bestpara=bestpara;
            out_struct.cv_full_result=result;
            out_struct.cv_best_result=min(min(result));
            varargout={bestpara,out_struct};
        end
    end
end