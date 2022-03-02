function [dcor,p]= partialdistcorr(x,y,numIter,verbose,pfor,nlin)
% [dcor,p]= distcorr(x,y,numIter,verbose,pfor)
%
% This function gets the relationship between the two columns in x
% (distance or pearson correlation) while controlling for all variables in
% y. x and y can either be tables or arrays. By "controlling" I mean that
% partialdistcorr will regress out variables in y and perform correlation
% over the residuals. 
%
% numIter : for distance correlation, a p-value is computed by permuting
% your data numIter times. Set to 0 to avoid doing permutations. I
% recommend setting to 5000. 
%
% verbose : if set to true, each permutation will be printed to display so
% you can track progress
%
% pfor : if set to true, permutations will be executed in parallel
%
% nlin : if set to true, distance correlation will be computed. If set to
% false, Pearson correlation will be computed (ommiting any rows w/NaNs). 
%
%
% Alex Teghipco // alex.teghipco@sc.edu

if ~istable(x)
    x = array2table(x);
end
if ~istable(y)
    y = array2table(y);
end
tab = [x y];

mdl1 = fitlm(tab,'interactions','ResponseVar',tab.Properties.VariableNames{1},...
    'PredictorVars',tab.Properties.VariableNames(3:end));
mdl2 = fitlm(tab,'interactions','ResponseVar',tab.Properties.VariableNames{2},...
    'PredictorVars',tab.Properties.VariableNames(3:end));

if nlin
    [dcor,p]= distcorr(mdl1.Residuals.Raw,mdl2.Residuals.Raw,numIter,verbose,pfor);
else
    [dcor,p]= corr(mdl1.Residuals.Raw,mdl2.Residuals.Raw,'rows','pairwise');
end