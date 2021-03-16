%% Demo

set(groot, 'defaultAxesTickLabelInterpreter','none');
set(groot,'defaulttextinterpreter','none');
set(groot, 'defaultLegendInterpreter','none');
set(groot, 'defaulttextfontsize',18);
set(groot, 'defaultLegendfontsize',18);
set(groot, 'defaultAxesfontsize',18);

%% Load dataSet info
dSID = 21190;

testC = DataSetDAN(dSID);

%% Waterfall plot 1

diag = 'IP2A';
sum1 = @(x) sum(x,2);

testC.waterfallPlot(diag, sum1);

%% Correlation plot 1

FSA2 = {'BPMS_LI20_3315_X'};
testC.correlationPlot(FSA2);


%% sort data

testC.sortOnFSArray(FSA2);
testC.correlationPlot(FSA2);


%% Waterfall plot 2

testC.waterfallPlot(diag, sum1);


%% select data
selFcn = @(x) (x > -0.3 & x < -0.2);

testC.selectImages(FSA2,selFcn);

testC.correlationPlot(FSA2);



%% Reset selection and sort

testC.selectImages;
testC.sortOnFSArray;

testC.correlationPlot(FSA2);

