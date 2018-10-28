function Ptable(tab,pP,times,save_name)

% make the table of p-value of test statistics 

% input:
% pP		-		p-values from WiCHProb function
% times		-		number of time, 4th dimension of the data
% 
% output:
% a table for different mean  P-value
if isempty(tab)
	tab = uifigure;
    tab.Name = 'P-values';
    tab.Position = [100,100,650,310];
end
saveP = round(pP,4);

colnames = strsplit(sprintf('t%d = ... = t%d\n',[1:times-2;repmat(times,[1,times-2])]),'\n');
colnames{times-1} = sprintf('t%d = t%d',times-1,times);
newP = num2cell(saveP,1);

str = strsplit(sprintf('t%d = t%d\n',[1:times-1;2:times]),'\n');
str = str([end, 1:end-1]);
rownames = str;
rownames{1} = 'Omnibus';

%tdata = table(newP{:},'VariableNames',colnames,'RowNames',rownames);
tdata = table(newP{:});

uit = uitable('Parent',tab,'Data',tdata,...
            'ColumnName',colnames,...
            'ColumnWidth','auto',...
            'RowName',rownames,...
            'Position',[10,10,600,300]);
        
%save the table into text file
colN = strsplit(sprintf('t%d_to_t%d\n',[1:times-2;repmat(times,[1,times-2])]),'\n');
colN{times-1} = sprintf('t%d_to_t%d',times-1,times);
str = strsplit(sprintf('t%d_equal_t%d\n',[1:times-1;2:times]),'\n');
str = str([end, 1:end-1]);
rowN = str;
rowN{1} = 'Omnibus';

newP1 = num2cell(saveP,1);
tdata1 = table(newP1{:},'VariableNames',colN,'RowNames',rowN);
name = [save_name,'.txt'];
writetable(tdata1,name,'WriteRowNames',true,'Delimiter',' ')  

end