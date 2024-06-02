clear;
clc;


addpath('data');
addpath('function');
fileList = dir('data\*.mat');
for i = 1:length(fileList)
    fileName = fileList(i).name;
    run(fileName);
end