
clc
clear all

%
A = csvread('results.csv');
plot(A([2:end],1),A([2:end],2),'linewidth',2)
xlabel('iterations')
ylabel('cost')