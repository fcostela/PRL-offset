%bceas = [];
%average = [];
% 1 NV 2 CVL (N=20) 3 Defocus (N=15)

% need to load them individually, because the var has the same name
load('nv_hist_allDegreesGOOD.mat')
load('defocus_hist_allDegrees.mat')
load('cvl_histDegreesEqual.mat')


j = 3;

for i=1:20
    
    values = find (~offset_mean{i}(:,1));
    if ~isempty(values)
        offset_mean{i}(values,1)= nan;
        offset_mean{i}(values,2)= nan;
    end
    x = offset_mean{i}(:,1);    
    y = offset_mean{i}(:,2);
    bceas(i,j) = bcea(x,y);
    c = sqrt(x.^2 + y.^2);   
    average(i,j) = nanmean(c);
end


% creating defocus levels matrix with offset
load('defocus_hist_allDegrees.mat');
load('./defocus/defocusdata.mat');
for i=1:15
        
    values = find (~offset_mean{i}(:,1));
    if ~isempty(values)
        offset_mean{i}(values,1)= nan;
        offset_mean{i}(values,2)= nan;
    end
    
    for j=1:length(offset_mean{i})        
        x = offset_mean{i}(j,1);    
        y = offset_mean{i}(j,2);     
       value(i,j) = sqrt(x.^2 + y.^2);   
    end
end

for i=1:5
    defocus{i} = [];
end
for i=1:15
        
%     values =find (~offset_mean{i}(:,1));
%     if ~isempty(values)
%         offset_mean{i}(values,1)= nan;
%         offset_mean{i}(values,2)= nan;
%     end     
  
    for j=1:5
        cols = find( defocuslevel(i,:)==j);      
        defocus{j} = [defocus{j};offset_mean{i}(cols,:)];
        
    end
end

for k = 1:5
    
    x = defocus{k}(:,1);    
    y = defocus{k}(:,2);     
    defocus_bcea(k) =  bcea(x,y);
end




k = 1;
for i=1:15
    for j=1:20
        level_offset(k,1) = value(i,j);
        level_offset(k,2) = defocuslevel(i,j);
      %  level_offset(k,3) = bceas(i,j);
        k = k+1;
    end
end
