
function[bcea1]=bcea(X,Y)
X=X';
p = 0.68;

chisq = -2*log(1-p); %for a probability of 68.2%
%for X coordinates of the fixation data
%Y= dictLookup(dict, 'Y');
Y=Y';
% Find number of samples in the X array.
n=size(X,2);
%%

% THIS SECTION CALCULATES THE PEARSON PRODUCT MOMENT CORRELATION, r, AND
% BIVARIATE CONTOUR ELLIPSE AREA (BCEA)
% Calculate means and standard deviations of X and Y
meanx=nanmean(X);
meany=nanmean(Y);
sdx=nanstd(X);
sdy=nanstd(Y);
n = size(X,2);
% Calculate sum(X =x), sum(y), (sum(x))^2=sxsq, and (sum(y))^2=sysq
sumx=nansum(X);
sumy=nansum(Y);
sxsq=sumx^2;
sysq=sumy^2;
% Calculate sum(xy)
xy=X.*Y;
sumxy=nansum(xy);
% Calculate sum x^2 and sum y^2
xsq=X.^2;
sumxsq=nansum(xsq);
ysq=Y.^2;
sumysq=nansum(ysq);
% Calculate r: Formula broken down into A, B, and C for convenience. r=A/(B*C)
A=(n*(sumxy))-((sumx)*(sumy));
B=sqrt((n*(sumxsq))-(sxsq));
C=sqrt ((n*(sumysq))-(sysq));
r=A/(B*C);
% Calculate BCEA
bcea1=2*pi*(chisq/2)*sdx*sdy*(sqrt(1-r^2));
%%
end