function [dcor,p]= distcorr(x,y,numIter,verbose,pfor)
% [dcor,p]= distcorr(x,y,numIter,verbose,pfor)
%
% This function calculates the distance correlation between x and y.
%
% Specify numIter as 0 in order to avoid permutation analysis, otherwise
% numIter determines the number of permutations that will be done to
% produce a p-value. 
%
% verbose can be set to a logical true or false to print each iteration as
% it is completed. 
%
% pfor can be set to a logical true or false to perform permutations in
% parallel.
%
% Reference: http://en.wikipedia.org/wiki/Distance_correlation 
%
% Date: 18 Jan, 2013
% Author: Shen Liu (shen.liu@hotmail.com.au).
%
% Date: 2 March, 2022
% Alex Teghipco (alex.teghipco@sc.edu) added
% permutation testing and parfor options.

% Check if the sizes of the inputs match
if size(x,1) ~= size(y,1);
    error('Inputs must have the same number of rows')
end

% Delete rows containing unobserved values
N = any([isnan(x) isnan(y)],2);
x(N,:) = [];
y(N,:) = [];

% Calculate doubly centered distance matrices for x and y
a = pdist2(x,x);
mcol = mean(a);
mrow = mean(a,2);
ajbar = ones(size(mrow))*mcol;
akbar = mrow*ones(size(mcol));
abar = mean(mean(a))*ones(size(a));
A = a - ajbar - akbar + abar;

b = pdist2(y,y);
mcol = mean(b);
mrow = mean(b,2);
bjbar = ones(size(mrow))*mcol;
bkbar = mrow*ones(size(mcol));
bbar = mean(mean(b))*ones(size(b));
B = b - bjbar - bkbar + bbar;

% Calculate squared sample distance covariance and variances
dcov = sum(sum(A.*B))/(size(mrow,1)^2);

dvarx = sum(sum(A.*A))/(size(mrow,1)^2);
dvary = sum(sum(B.*B))/(size(mrow,1)^2);

% Calculate the distance correlation
dcor = sqrt(dcov/sqrt(dvarx*dvary));

% Do a permutation test for p-value...added by Alex
mrowy=mrow;
if numIter > 0
    rng shuffle
    if ~pfor
        for i = 1:numIter
            if verbose
                disp(['Working on p-value calcuation...' num2str(i) ' of ' num2str(numIter)]);
            end
            xP = x(randperm(length(x)));
            
            a = pdist2(xP,xP);
            mcol = mean(a);
            mrow = mean(a,2);
            ajbar = ones(size(mrow))*mcol;
            akbar = mrow*ones(size(mcol));
            abar = mean(mean(a))*ones(size(a));
            A = a - ajbar - akbar + abar;
            
            % Calculate squared sample distance covariance and variances
            dcov = sum(sum(A.*B))/(size(mrowy,1)^2);
            
            dvarx = sum(sum(A.*A))/(size(mrowy,1)^2);
            %dvary = sum(sum(B.*B))/(size(mrowy,1)^2);
            
            % Calculate the distance correlation
            permD(i,1) = sqrt(dcov/sqrt(dvarx*dvary));
        end
    else
        parfor i = 1:numIter
            xP = x(randperm(length(x)));
            
            a = pdist2(xP,xP);
            mcol = mean(a);
            mrow = mean(a,2);
            ajbar = ones(size(mrow))*mcol;
            akbar = mrow*ones(size(mcol));
            abar = mean(mean(a))*ones(size(a));
            A = a - ajbar - akbar + abar;
            
            % Calculate squared sample distance covariance and variances
            dcov = sum(sum(A.*B))/(size(mrowy,1)^2);
            
            dvarx = sum(sum(A.*A))/(size(mrowy,1)^2);
            %dvary = sum(sum(B.*B))/(size(mrowy,1)^2);
            
            % Calculate the distance correlation
            permD(i,1) = sqrt(dcov/sqrt(dvarx*dvary));
        end
    end
    %p = 1-(length(find(dcor >= permD))/numIter);
    p = length(find(permD > dcor))/numIter;
else
    p = NaN;
end
