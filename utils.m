classdef utils
    methods(Static)
        %%   ��һ����ά������ҳ���Сֵ��λ�ã��������β��ȶ������ڵ���Сֵ
        function [row_id,col_id]=from_table(table,varargin)
            %   varargin�е�Ԫ�ؿ�����Сֵ�����ѡȡ��Ԫ������Ϊ��
            %       valid_timesָ���˶��ٱ�����Сֵ�ĵ�ᱻ��Ϊ�����ε㣬Ĭ��inf��һ����ʹ�����ι��ܣ���ȡ3
            %       horizontal_range��ʾ���ε���Χ������ٷ�ΧҲ��������������Ĭ��0�����ݺ���Χ����Ӧ���ѡȡ
            %       vertical_rangeͬ�ϣ�����Χ
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
            invalid_map=table>threshold;                        %   һ��"�Ƿ�ֵ"��flagͼ
            [invalid_row,invalid_col]=find(invalid_map);    %   �ҵ����е�"�Ƿ���"
            filtered_table=table;
            for i=1:length(invalid_row)
                %   ��ÿһ���Ƿ�����Χ�ó�inf
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

        %% ���ƶԽ���ͼ�����ڹ۲�Ԥ��Ч��(����model.evaluate)
        %   ����Ĳ��������Ƿ���ɢ���Ա߱�ע���
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
            
            %   д���
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
        
        %% ���Ƶ���������ͼ(����model.tune)
        function single_variable_min(data,x_range,x_label)
            figure;hold on
            plot(x_range,data);                     %   ������������
            
%             idx=find(data<1.05*min(min(data))); %   ���ƽӽ���Сֵ�ĵ�
%             scatter(x_range(idx),data(idx),'go')
            
            [minval,index]=min(data);
            scatter(x_range(index),minval,'ro')      %   Ѱ����Сֵ��
            xlabel(strrep(x_label,'_','\_'))
            ylabel('RMSE')
        end
        
        %%  ����˫��������ͼ(����model.tune)
        function double_variable_min(datamat,x_range,y_range,x_label,y_label)
            figure;hold on
            minvalue=min(min(datamat));     %   �Ƴ���ֵ���ߵĲ���
            datamat(datamat>3*minvalue)=nan;
            %%  ���Ƶȸ���ͼ
            [X,Y]=meshgrid(x_range,y_range);
            datamat=datamat';
            contourf(X,Y,datamat)
            %%  ����Сֵ
%             [y,x]=find(datamat<1.05*min(min(datamat)));   %   �����д���������1.05�����RMSE�ڵĽ⣬��ʾΪ�̵�
%             scatter(x_range(x),y_range(y),'go')

            [y,x]=find(datamat==min(min(datamat)));
            scatter(x_range(x),y_range(y),'ro')
            xlabel(strrep(x_label,'_','\_'))
            ylabel(strrep(y_label,'_','\_'))
        end
        
        %% �������-��������ͼ(����model.tune)
        function class_variable_min(datamat,classes,v_range,classes_name,v_label)
            figure
            hold on
            for i=1:length(classes)
                plot(v_range,datamat(i,:))
            end
            legend(classes_name)
            xlabel(strrep(v_label,'_','\_'))
        end
        
        %% ���������״����ͼ
        function class_min(datamat,classes)
            bar(datamat)
            set(gca, 'XTickLabel', classes)
            ylabel('RMSE')
        end
        
        %%  ����ά�Ȳ��ȵľ�����ӣ�cross validationʱ�п��������������⣩
        function result=add_mismatch(data1,data2)
            [rows1,cols1]=size(data1);
            [rows2,cols2]=size(data2);
            rows=min(rows1,rows2);
            cols=min(cols1,cols2);
            result=data1(1:rows,1:cols)+data2(1:rows,1:cols);
        end
        
        %% ����������ĵȸ���ͼ
        %   ������Ҫ��scatteredInterpolant
        function surf_scatter(x,y,z)
            intervals=20;   %   ������Ϊ���ɸ�������в�ֵ
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