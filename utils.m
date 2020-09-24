classdef utils
    methods(Static)
        %%   从一个二维表格中找出最小值的位置，可以屏蔽不稳定区域内的最小值
        function [row_id,col_id]=from_table(table,varargin)
            %   varargin中的元素控制最小值点如何选取，元素依次为：
            %       valid_times指定了多少倍于最小值的点会被认为是屏蔽点，默认inf；一般若使用屏蔽功能，可取3
            %       horizontal_range表示屏蔽点周围横向多少范围也会纳入屏蔽区域，默认0；根据横向范围与响应情况选取
            %       vertical_range同上，纵向范围
            [rows,cols]=size(table);
            if length(varargin)>0
                valid_times=varargin{1};
            else
                valid_times=inf;
            end
            if length(varargin)>1
                horizontal_range=varargin{2};
            else
                horizontal_range=0;
            end
            if length(varargin)>2
                vertical_range=varargin{3};
            else
                vertical_range=0;
            end
            min_value=min(min(table));
            threshold=min_value*valid_times;
            invalid_map=table>threshold;                        %   一个"非法值"的flag图
            [invalid_row,invalid_col]=find(invalid_map);    %   找到所有的"非法点"
            filtered_table=table;
            for i=1:length(invalid_row)
                %   把每一个非法点周围置成inf
                row=invalid_row(i);
                col=invalid_col(i);
                row_lb=max(1,row-vertical_range);
                row_ub=min(rows,row+vertical_range);
                col_lb=max(1,col-horizontal_range);
                col_ub=min(cols,col+horizontal_range);
                filtered_table(row_lb:row_ub,col_lb:col_ub)=inf;
            end
            if sum(sum(isinf(filtered_table)))==cols*rows
                error('Range too large, no valid value!')
            end
            [row_ids,col_ids]=find(filtered_table==min(min(filtered_table)));
            row_id=row_ids(1);
            col_id=col_ids(1);
        end

        %% 绘制对角线图，用于观察预测效果(用于model.evaluate)
        %   额外的参数控制是否在散点旁边标注序号
        function comparision(y_true,y_pred,varargin)
            figure;hold on
            scatter(y_true,y_pred)
            
            if length(varargin)>0
                label_flag=varargin{1};
            else
                label_flag=false;
            end
            
            maxval=max(max(y_true),max(y_pred));
            minval=min(min(y_true),min(y_pred));
            d=maxval-minval;
            axis([minval-0.1*d,maxval+0.1*d,minval-0.1*d,maxval+0.1*d])
            line([minval-0.1*d,maxval+0.1*d],[minval-0.1*d,maxval+0.1*d],'color','k','linestyle','--')
            
            %   写编号
            if label_flag
                for i=1:length(y_true)
                    text(y_true(i),y_pred(i),num2str(i),'HorizontalAlignment','Center')
                end
            end
            
            err=y_true-y_pred;
            rmse=sqrt(sum(sum(err.^2))/length(err));
            title(['RMSE = ',num2str(rmse)])
            xlabel('True','FontSize',9)
            ylabel('Pred','FontSize',9)
        end
        
        %% 绘制单变量调试图(用于model.tune)
        function single_variable_min(data,x_range,x_label)
            figure;hold on
            plot(x_range,data);                     %   绘制性能曲线
            
%             idx=find(data<1.05*min(min(data))); %   绘制接近最小值的点
%             scatter(x_range(idx),data(idx),'go')
            
            [minval,index]=min(data);
            scatter(x_range(index),minval,'ro')      %   寻找最小值点
            xlabel(strrep(x_label,'_','\_'))
            ylabel('RMSE')
        end
        
        %%  绘制双变量调试图(用于model.tune)
        function double_variable_min(datamat,x_range,y_range,x_label,y_label)
            figure;hold on
            minvalue=min(min(datamat));     %   移除掉值过高的部分
            datamat(datamat>3*minvalue)=nan;
            %%  绘制等高线图
            [X,Y]=meshgrid(x_range,y_range);
            datamat=datamat';
            contourf(X,Y,datamat)
            %%  找最小值
%             [y,x]=find(datamat<1.05*min(min(datamat)));   %   这两行代码用于找1.05倍最佳RMSE内的解，显示为绿点
%             scatter(x_range(x),y_range(y),'go')

            [y,x]=find(datamat==min(min(datamat)));
            scatter(x_range(x),y_range(y),'ro')
            xlabel(strrep(x_label,'_','\_'))
            ylabel(strrep(y_label,'_','\_'))
        end
        
        %% 绘制类别-变量调试图(用于model.tune)
        function class_variable_min(datamat,classes,v_range,classes_name,v_label)
            figure
            hold on
            for i=1:length(classes)
                plot(v_range,datamat(i,:))
            end
            legend(classes_name)
            xlabel(strrep(v_label,'_','\_'))
        end
        
        %% 绘制类别柱状调试图
        function class_min(datamat,classes)
            bar(datamat)
            set(gca, 'XTickLabel', classes)
            ylabel('RMSE')
        end
        
        %%  用于维度不等的矩阵相加（cross validation时有可能遇到这种问题）
        function result=add_mismatch(data1,data2)
            [rows1,cols1]=size(data1);
            [rows2,cols2]=size(data2);
            rows=min(rows1,rows2);
            cols=min(cols1,cols2);
            result=data1(1:rows,1:cols)+data2(1:rows,1:cols);
        end
        
        %% 不规则网格的等高线图
        %   可能需要用scatteredInterpolant
        function surf_scatter(x,y,z)
            intervals=20;   %   将划分为若干个网格进行插值
            xmin=min(x);xmax=max(x);
            ymin=min(y);ymax=max(y);
            xq=linspace(xmin,xmax,intervals);
            yq=linspace(ymin,ymax,intervals);
            z3 = griddata(x,y,z,xq,yq,'natural');
            figure
            plot3(x,y,z,'mo')
            hold on
            mesh(xq,yq,z3)
        end

    end
end