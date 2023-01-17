function [bcpR, bcpP, pR, pP, pRM, pPM] = pdc(x,y,z,corrType,incMat)
% This function implements the recursive formula to compute *any* nth order
% two-sample partial correlation or partial bias corrected distance
% correlation. See Szekely's lecture
% slides for more information/justification for this approach:
% https://stat.wisc.edu/wp-content/uploads/sites/870/2020/03/SzekelyGabor.pdf
%
% x must be an n x 1 vector that you are trying to relate to y (also an n x
% 1 vector), where n is the number of samples in your data. z is an n x p
% matrix with p nuissance variables that you would like to condition on.
% corrType can be 'pearson' or 'distance'. incMat can be set to true to
% compute partial correlation coefficient/p-value using matlab's internal
% stats toolbox function (i.e., without using the recursive formula).
%
% NOTE: matlab's partial r coefficient will deviate slightly at higher
% order partial correlations due to precision (i.e., 4th order and above).
%
% Outputs: 
% bcpR: bias corrected partial distance correlation coefficient
% bcpP: bias corrected partial distance correlation p-value
% pR: partial correlation coefficient
% pP: partial correlation p-value
% pRM: partial correlation coefficient from matlab's internal partialcorr
% (i.e., not using recursive formula; uses regression method)
% pPM: partial correlation p-value from matlab's internal partialcorr
%
% [bcpR, bcpP, pRM, pCM,] = pdc(x,y,z,corrType)
% Example call for partial distance correlation: [bcpR, bcpP, pRM, pCM,] = pdc(x,y,z,'distance',false)
%
% alex.teghipco@sc.edu 

%% startup
if isempty(incMat)
    incMat = true;
    disp('No incMat argument...setting to true')
end
if isempty(corrType)
    corrType = 'pearson';
     disp('No corrType argument...setting to pearson')
end
bcpR = [];bcpP = [];pR = [];pP = [];pRM = [];pPM = [];
if isempty(which('bcdistcorr'))
   tmp = which('pdc.m');
   addpath(genpath(tmp)) ;
   disp('You do not have bcdistcorr.m in your path...attempting to add it by assuming it is in the same directory as this function...')
end

d = [x y z];
n = size(d,2)-2; % order eg 2nd

%% Sets up cell structure and gets correct indices to correlate
lvl = n;
while lvl > 0
    if lvl == n
        idx = [1:2 3:3+(lvl-1)];
        [tmpl] = getord([idx(3:end)]);
        mz{lvl}{1}(1,:) = idx(1:end-1);
        mz{lvl}{1}(2,:) = [idx(1) tmpl];
        mz{lvl}{1}(3,:) = [idx(2) tmpl];
    else
        tmp = vertcat(mz{lvl+1}{:});
        for j = 1:size(tmp,1)%1:sum(cellfun('size',mz{lvl+1},2))%1:size(mz{lvl+1},1)
            idx = tmp(j,:);%mz{lvl+1}{1}(j,:);
            if lvl == 1
                [tmpl] = getord(idx(3:end));
            else
                [tmpl] = getord(idx(3:end));
            end
            mz{lvl}{j}(1,:) = idx(1:end-1);
            mz{lvl}{j}(2,:) = [idx(1) tmpl];
            mz{lvl}{j}(3,:) = [idx(2) tmpl];
        end
    end
    lvl = lvl - 1;
end

%% now do the correlations
for i = 1:n
    clear tmpr
    tmp = vertcat(mz{i}{:});
    if i == 1
        for j = 1:size(tmp,1)
            if strcmpi(corrType,'pearson')
                tmpr(j,1)= corr(d(:,tmp(j,1)),d(:,tmp(j,2)));
            elseif strcmpi(corrType,'distance')
                [tmpr(j,1), ~, ~, ~] = bcdistcorr(d(:,tmp(j,1)),d(:,tmp(j,2)));
            end
        end
        tmp2 = reshape(tmpr,size(tmp,1)/size(mz{i},2),size(mz{i},2))';
    else
        for j = 1:size(tmp,1)
            tmpr(j,1) = (rk{i-1}{j}(1) - (rk{i-1}{j}(2)*rk{i-1}{j}(3)))/sqrt(((1-(rk{i-1}{j}(2)).^2))*((1-(rk{i-1}{j}(3)).^2)));
        end
        tmp2 = reshape(tmpr,size(mz{i},2),size(tmp,1)/size(mz{i},2));
    end
    for j = 1:size(tmp2,1)
        rk{i}{j} = tmp2(j,:);
    end
end

% hack for now--in case 1st order try, in case > 1 catch...
try
    o = (rk{end}{1}(1) - (rk{end}{1}(2))*(rk{end}{1}(3)))/sqrt(((1-(rk{end}{1}(2)).^2))*((1-(rk{end}{1}(3)).^2)));
catch
    o = (rk{end}{1} - (rk{end}{2})*(rk{end}{3}))/sqrt(((1-(rk{end}{2}).^2))*((1-(rk{end}{3}).^2)));
end
if o > 1
    o = 1;
end
if strcmpi(corrType,'pearson')
    pR = o;
    df = max(size(x,1)-size(z,2)-2,0);
    t = sign(pR) .* Inf;
    k = (abs(pR) < 1);
    t(k) = pR(k) ./ sqrt(1-pR(k).^2);
    t = sqrt(df).*t;
    pP = 2*tcdf(-abs(t),df);
elseif strcmpi(corrType,'distance')
    bcpR = o;
    n = size(x,1);
    %M = n*(n-3)/2; % for normal distance correlation 
    M = n*(n-4)/2; % for parital distance
    T = sqrt(M-1) * bcpR / sqrt(1-bcpR^2);
    df = M-1;
    bcpP = 1 - tcdf(T, df);
end
if incMat
    [pRM,pPM] = partialcorr(x,y,z);
end
