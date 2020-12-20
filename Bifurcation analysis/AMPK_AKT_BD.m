file_name = '/path/to/parameter/set/file';

data = readtable(file_name);
par_sets = table2cell(data);
par_set = cell2mat(par_sets(:,1:25));

for i = 1:size(par_set,1)
    parameters = par_set(i,:);
    [x,v,s,h,f] = AMPK_AKT_Bifurcation(parameters, 5); 
    a = x(5,:); %bifurcation parameter
    b = x(2,:); %AKT
    c = x(1,:); %AMPK
    plot(a,c,'b-','LineWidth',2);
    xlabel('kac_AMPK');
    ylabel('AMPK');
    title('Bifurcation of AMPK levels ');
    xlim([0 2])
    ylim([0 100])
    y=x';
    %save('Stat1_Stat3_miR200_bifurcation.txt','y','-ascii');

    dat=struct2table(s);
    % %% Based on eigenvalues to judge stable vs. unstable states
    ind = table2array(dat(:,1));


    %%
    fig= figure('Color',[1 1 1],'units','normalized');%,'outerposition',[0 0 1 1]);
    grid on
    for i=1:length(ind)-1
        if rem(i,2)
            plot(a(ind(i):ind(i+1)),c(ind(i):ind(i+1)),'-bs','LineWidth',2,'MarkerIndices',[ind(1),ind(i)]);
            hold on;
        else
            plot(a(ind(i)+1:ind(i+1)),c(ind(i)+1:ind(i+1)),'-rs','LineWidth',2,'MarkerIndices',[ind(1),ind(i)]);
            hold on;
        end
    end
    xlim([0 2]);
    ylim([0 100]);
    xlabel('k_{ac} AMPK');
    ylabel('AMPK');
    title('Bifurcation of AMPK levels ');
    hold on;
    % saveas(fig,'Bifurcation of STAT3 levels','jpeg')
    % saveas(fig,'Bifurcation of STAT3 levels','pdf')
    % savefig('Bifurcation of STAT3 levels')
end
