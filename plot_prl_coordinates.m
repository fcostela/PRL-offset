% load the offset knowns and plot against the PRL coordinates

%load prl_coords
subjects_vf = {'1208dd', '1732ri', '2484ao', '3024mn', '1783rz' , '2765ma', '0212dy', '2338ms'};
subjects_gf = {'2765ma', '0212dy'};%, '2338ms'};

% This code was helpful to convert Goldman data points from a text file
% with 2 columns into matlab

fid = fopen('0212dyOD_Goldman.txt');
D = textscan(fid, '%s','delimiter', '\n');
data = [];
for i=2:length(D{1})
     if ~isempty(sscanf(cell2mat(D{1}(i)), '%g, %g ',1))
            data(i-1,:) = sscanf(cell2mat(D{1}(i)), '%g, %g ',2)';
     end
     
end

offset_vf_x = [-9.791, -7.862714748, -6.255101563, -6.179855382, -1.388185969  , 1.731791212, -0.92229183, -2.29716367];
offset_vf_y = [-7.954, -7.843919581, -3.980085017, -4.51467952, 9.374760004 , 1.214092356, 1.565788548, 8.627981295];
offset_gf_x = [1.731791212, -0.92229183]; %, -2.29716367];
offset_gf_y = [1.214092356, 1.565788548];%, 8.627981295];

size = 150;
sizeScot = 50;
limits = 20;
for k=1:length(subjects_vf)
    j = 1;
    s = subjects_vf{k}
    load(s);
    
     subplot(4,2,k);
     
    if k<6
        valid = [];
       
%         for i=1:length(data)
%             if sqrt(data(i,1)^2 + data(i,2)^2) < 200
%                 valid(j) = i;
%                 j = j+1;
%             end
%         end
        
       % points = data(valid,:);
        points = data;
        
        datapoints =pix2deg(points,screenDist,screenSize,screenRes);
    else
        datapoints = data;
    end
    
    %transform pixels to degrees
    % screenDist=45;         distance in cm
    % screenSize=[60 34];    [horz vert] screen dimensions in cm
    % screenRes=[2560 1440]; screen dims in pixels    
    colortype = [0 0 0] ;
    colortype2 = [0 0 1];
    colortype3 = [1 0 0] 
    if k <8
    scatter(datapoints(:,1), datapoints(:,2), sizeScot,[colortype3(1) colortype3(2) colortype3(3)],'h' , 'filled');
    end
   % fill(datapoints(:,1), datapoints(:,2), [1 0 0]);
    hold on;
    
    plot(0,0,'+');
    
  
    scatter(offset_vf_x(k), offset_vf_y(k), size, [colortype2(1) colortype2(2) colortype2(3)],'o','filled');
      scatter(offset_vf_x(k), offset_vf_y(k), size+300, [colortype3(2) colortype3(2) colortype3(3)],'+');
    
    xlabel('Horizontal  (degrees)','FontName', 'Arial', 'FontSize', 18);
    ylabel('Vertical (degrees)','FontName', 'Arial', 'FontSize', 18);
    if k<6
        typeR = ' from Visual Fields';
    else
        typeR = ' from Goldman';
    end
    
    title([subjects_vf{k} typeR], 'FontSize', 20);
    line([-limits limits], [0 0], 'Color', 'k');
    line([0 0], [-limits limits], 'Color', 'k');    
    xlim([-limits limits]);ylim([-limits limits]);
    set(gca,'FontName', 'Arial', 'FontSize', 16);
   
end

% j =1;
% for i=1:2:length(coords)
%    subplot(7,2,j);
%    size = 100;
%    offset_x = coords(i,1);
%    offset_y = coords(i,2);
%    
%    prl_r_x = coords(i,3);
%    prl_r_y = coords(i,4);
%    
%    prl_l_x = coords(i+1,3);
%    prl_l_y = coords(i+1,4);
%    
%    va_l = coords(i+1,5);
%    va_r = coords(i,5);
% 
%    mymin_x = min([ offset_x prl_r_x prl_l_x]);
%    mymin_y = min([ offset_y prl_r_y prl_l_y]);
%    mymax_x = max([ offset_x prl_r_x prl_l_x]);
%    mymax_y = max([ offset_y prl_r_y prl_l_y]);
%    
%    scatter(offset_x, offset_y, size,[colortype(1) colortype(2) colortype(3)],'+');
%    hold on;
%    if va_r      
%        size = 300;
%    else
%        size = 100;
%    end
%    scatter(prl_r_x, prl_r_y, size, [colortype2(1) colortype2(2) colortype2(3)],'o','filled');
%    if va_l
%        size = 300;
%    else
%        size = 100;
%    end
%    scatter(prl_l_x, prl_l_y, size, [colortype3(1) colortype3(2) colortype3(3)],'o','filled');
%    j = j +1;       
%     xlabel('Horizontal  (degrees)','FontName', 'Arial', 'FontSize', 18);
%  ylabel('Vertical (degrees)','FontName', 'Arial', 'FontSize', 18);
% 
% line([-30 25], [0 0], 'Color', 'k');
% line([0 0], [-25 25], 'Color', 'k');
% 
% set(gca,'xlim',[mymin_x-5 mymax_x+5]);
% set(gca,'ylim',[mymin_y-5 mymax_y+5], 'FontName', 'Arial', 'FontSize', 18);
% set(gca, 'FontName','Arial','FontSize',16);
% end


