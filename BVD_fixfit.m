%THIS PROGRAM CALCULATES BIVARIATE ELLIPSES FOR NORMAL SUBJECTS' FOVEAL
%FIXATION DATA IN DEGREES^2.  THE INPUT IS UN-DISTORTED X,Y FOVEAL POSITION
%RELATIVE TO THE OPTIC DISC IN DEGREES FROM EXCEL SPREADHEET 'FIXATION RE-MEASURE'
%the BVD program modified for fixation fitting by AR-05/09
function[fovea]=BVD_fixfit(imageDir,nidekFileStructs,m,fovea,ecc_Theta,ecc_R)
for ctr = 1:length(nidekFileStructs)

fullFilename = [imageDir '\' nidekFileStructs(ctr).name];
% p=0.682;
p=0.95;
Response = questdlg('The default probability value for BVD is 0.682. Do you wish to change to 0.90?','Skip','Yes','No','Cancel','Yes');
if strcmp(Response, 'No')
    disp('68.2%-Continue')
end
if strcmp(Response, 'Yes')
    disp('90%');
    p=0.9;
end
if strcmp(Response, 'Cancel')
    disp('Quitting on User Request');

end

p = 0.68;

chisq = -2*log(1-p); %for a probability of 68.2%

% chisq = 2.279; %for a probability of 68.2%

    textFilename = convertImageFilenameToNidekFilename(fullFilename);
    [dict, pixperdeg] = parseLITextFile(textFilename);
    
%%
% Convert the points, since they appear in degrees and we want pixels.

X = dictLookup(dict, 'X');
%converting the column vector to row vector
X=X';

%for X coordinates of the fixation data
Y= dictLookup(dict, 'Y');
Y=Y';
% Find number of samples in the X array.
n=size(X,2);
%%

% THIS SECTION CALCULATES THE PEARSON PRODUCT MOMENT CORRELATION, r, AND
% BIVARIATE CONTOUR ELLIPSE AREA (BCEA)
% Calculate means and standard deviations of X and Y
meanx=mean(X);
meany=mean(Y);
sdx=std(X);
sdy=std(Y);
n = size(X,2);
% Calculate sum(X =x), sum(y), (sum(x))^2=sxsq, and (sum(y))^2=sysq
sumx=sum(X);
sumy=sum(Y);
sxsq=sumx^2;
sysq=sumy^2;
% Calculate sum(xy)
xy=X.*Y;
sumxy=sum(xy);
% Calculate sum x^2 and sum y^2
xsq=X.^2;
sumxsq=sum(xsq);
ysq=Y.^2;
sumysq=sum(ysq);
% Calculate r: Formula broken down into A, B, and C for convenience. r=A/(B*C)
A=(n*(sumxy))-((sumx)*(sumy));
B=sqrt((n*(sumxsq))-(sxsq));
C=sqrt ((n*(sumysq))-(sysq));
r=A/(B*C);
% Calculate BCEA
bcea1=2*pi*(chisq/2)*sdx*sdy*(sqrt(1-r^2));
%%

% THIS SECTION CALCULATES THE ELLIPSE GEOMETRY
% Calculate angles of ellipse axes. major axis slope=slpmaj; minor axis
% slope (slopmin)
beta=0.5*(atan((2*r*sdx*sdy)/((sdx^2)-(sdy^2))));
slpmaj=tan(beta);
slpmin=-(1/slpmaj);
%Calculate intersection of major and minor axes and ellipse.
% Minor axis: x1 (xone) and y1 (yone)
D=chisq*(1-(r^2))*((sdx*sdy)^2);
E=((slpmin^2)*(sdx^2))+(sdy^2)-(2*r*sdx*sdy*slpmin);
xone=meanx+(sqrt((D/E)));
yone=meany+(slpmin*(xone-meanx));
% Major axis: x2 (xtwo) and y2 (ytwo)
G=((slpmaj^2)*(sdx^2))+(sdy^2)-(2*r*sdx*sdy*slpmaj);
xtwo=meanx+(realsqrt((D/G)));
ytwo=meany+(slpmaj*(xtwo-meanx));

min_axis_inter_1=[xtwo ytwo]; %minor axis intersecting the ellipse
maj_axis_inter_1=[xone yone];%major axis intersecting the ellipse

%Calculate the semi-lengths of the major and minor axes. Lone (L1) is the
%length of the semi-minor ellipse axis, and L2 (Ltwo) is the semi-major
%length.
Lone=sqrt(((xone-meanx)^2)+((yone-meany)^2));
Ltwo=sqrt(((xtwo-meanx)^2)+((ytwo-meany)^2));
bcea2=pi*Lone*Ltwo;
Lmin=Lone*2;
Lmaj=Ltwo*2;


%finding the equation of the major axis
%(y-y1)=m(x-x1)
% ymaj-yone=slpmaj*(xmaj-xone);
constant_maj=-slpmaj*xone+yone;
% ymaj=(slpmaj*xmaj)+constant_maj;

%finding the equation of the minor axis
% ymin-ytwo=slpmin*(xmin-xtwo);
constant_minor=-slpmin*xtwo+ytwo;
%ymin=(slpmin*xmin)+constant_minor;
%%

%getting the coordinates of the major and minor axis intersection with the BVD

if (xone>meanx)
    xthree=meanx-(abs(meanx-xone));
else
    xthree=meanx+(abs(meanx-xone));
end
if (yone>meany)
    ythree=meany-(abs(meany-yone));
else
    ythree=meany+(abs(meany-yone));
end
if (xtwo>meanx)
    xfour=meanx-(abs(meanx-xtwo));
else
    xfour=meanx+(abs(meanx-xtwo));
end
if (ytwo>meany)
    yfour=meany-(abs(meany-ytwo));
else
    yfour=meany+(abs(meany-ytwo));
end

maj_axis_inter_2=[xthree ythree]; %major axis intersecting the ellipse
min_axis_inter_2=[xfour yfour];%minor axis intersecting the ellipse
center_BVD=[meanx meany];
Cof_maj_axis=[slpmaj,1];
Cof_min_axis=[slpmin,1];

%getting back assumed center of the fovea in temporal values
fovea(1)=(fovea(1)-1)/pixperdeg;
fovea(2)=(fovea(2)-1)/pixperdeg;
fovea=[fovea(1) fovea(2)];
%theta and eccentricity of the vector joining the center of the BVD and the
% assumed center of the fovea
[eccTheta_1 eccR_1] = cart2pol(meanx-fovea(1),meany-fovea(2)); 
% [eccTheta_1 eccR_1] = cart2pol(fovea(1)-meanx,fovea(2)-meany); 

%finding the angle between the major axis of BVD and the Fovea-PRL vector
beta=eccTheta_1-slpmaj;
%%

% Display ellipse parameters into a text file
outputDir = imageDir;
output_Filename= strrep(textFilename,'.mfd','.txt'); 
output_Filename_test= strrep(textFilename,'.mfd','_test.txt'); 
%   outputFilename  = [outputDir '\' output_Filename_1];
  output_Imagename=strrep(textFilename,'.mfd','.emf');
%   outputImagename = [outputDir '\' textFilename];
  
test_datafile(output_Filename_test,p,chisq,center_BVD,min_axis_inter_1,...
    min_axis_inter_2,maj_axis_inter_1,maj_axis_inter_2,slpmaj, sdx,sdy,...
    r,bcea1,bcea2,Lmaj,Lmin,Cof_maj_axis,Cof_min_axis,...
    fovea,ecc_Theta,ecc_R,eccTheta_1,eccR_1,beta);

s=struct('p',p,'center_BVD',center_BVD,'maj_axis_inter_1',maj_axis_inter_1,...
         'maj_axis_inter_2',maj_axis_inter_2,'min_axis_inter_1',...
          min_axis_inter_1,'min_axis_inter_2',min_axis_inter_2,'slpmaj',slpmaj,'sdx',sdx,...
         'sdy',sdy,'bcea1',bcea1,'bcea2',bcea2,'Lmaj',Lmaj,'Lmin',Lmin,...
         'Cof_maj_axis',Cof_maj_axis,'Cof_min_axis',Cof_min_axis,'fovea',fovea,...
         'ecc_Theta',ecc_Theta,'ecc_R',ecc_R,'eccTheta_1',eccTheta_1,'eccR_1',eccR_1,...
         'beta',beta);
     
SPSSvariableStore_AR(imageDir,output_Filename,'initialize',s);
SPSSvariableStore_AR(imageDir,output_Filename,'getInstance',s);
SPSSvariableStore_AR(imageDir,output_Filename,'add',s);
SPSSvariableStore_AR(imageDir,output_Filename,'output',s);
%% Plotting

%converting the points to pixels
X = round(X * pixperdeg)+1;
Y = round(Y * pixperdeg)+1;

%plotting the X and Y fixation points
figure(m);
hold on
for i=1:n
    plot (X(i),Y(i),'or','MarkerFaceColor',[.49 1 .63],'MarkerSize',6);
end

hold on
meanx_1=(meanx * pixperdeg)%+1;
meany_1=(meany * pixperdeg)%+1;
xone=(xone * pixperdeg)%+1;
yone=(yone * pixperdeg)%+1;
xtwo=(xtwo * pixperdeg)%+1; 
ytwo=(ytwo * pixperdeg)%+1;
plot (meanx_1,meany_1,'+c','MarkerSize',12);

% Plot minor axis
if (xone>meanx_1)
    xthree=meanx_1-(abs(meanx_1-xone));
else
    xthree=meanx_1+(abs(meanx_1-xone));
end
if (yone>meany_1)
    ythree=meany_1-(abs(meany_1-yone));
else
    ythree=meany_1+(abs(meany_1-yone));
end
a=[xone xthree];
b=[yone ythree];

hold on

% Plot major axis
if (xtwo>meanx_1)
    xfour=meanx_1-(abs(meanx_1-xtwo));
else
    xfour=meanx_1+(abs(meanx_1-xtwo));
    xfour=(xfour * pixperdeg)%+1;
end
if (ytwo>meany_1)
    yfour=meany_1-(abs(meany_1-ytwo));
else
    yfour=meany_1+(abs(meany_1-ytwo));
end
plot (xone,yone,'oy','LineWidth',8)
plot (xthree,ythree,'og','LineWidth',8)
plot (xtwo,ytwo,'om','LineWidth',8)
plot (xfour,yfour,'ob','LineWidth',8)
plot (a,b,'m','LineWidth',2);
c=[xtwo xfour];
d=[ytwo yfour];
hold on
plot (c,d,'m','LineWidth',2);
%%

% Plott BVD using 'Ellipse' program

%  Ltwo=(Ltwo * pixperdeg)+1;
%  Lone=(Lone * pixperdeg)+1;
%  slpmaj=(slpmaj * pixperdeg)+1;
% ellipse_AR shifted the ellipse for some reason...
%ellipse_AR(Ltwo,Lone,slpmaj,meanx,meany,'m',pixperdeg,500)
ellipse_Plot(Ltwo,Lone,slpmaj,meanx,meany,'m',pixperdeg,500)


axis equal %forces to use the same scale on the X and Y axes
xlabel('Fixation Data (X coordinate)');
ylabel('Fixation Data (Y Coordinate)');
saveas(gca,output_Imagename); %saving the output plot as emf file

end
end
