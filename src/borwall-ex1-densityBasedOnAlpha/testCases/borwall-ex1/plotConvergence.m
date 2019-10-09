
clc
clear all

%
%A = csvread('results.csv');
%plot(A([2:end],1),A([2:end],2))
%xlabel('iterations')
%ylabel('cost')

alphaAbsMin = 1e-8;
q = 0.01;
alphaAbsMax = 0.5;


rho  = @(alpha)(q*(alpha - alphaAbsMax)./( (alphaAbsMin - alphaAbsMax)*(1 + q) - (alpha - alphaAbsMax)));

a = linspace(alphaAbsMin,alphaAbsMax,100);
plot(a,rho(a))