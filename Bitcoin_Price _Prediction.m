bitcoin_data = readtable('NSE_data.csv');
open = bitcoin_data.Open;
close = bitcoin_data.Close;
high= bitcoin_data.High;
adj_close=bitcoin_data.AdjClose;
low=bitcoin_data.Low;
date=datetime(bitcoin_data.Date);
vol=bitcoin_data.Volume;
f_start=0.5;
r=open;
n=length(r);
n1=fix(n*f_start);  %estimation sample
n2=n-n1;

%For h = 1
%aggregate dependent variable (containing true values)
h = 1;
y=zeros(n-h+1,1);
for i=1:(n-h+1);
    s=0;
    for j=0:(h-1);s=s+r(i+j);end;
    y(i)=s;
end;
y_true=y((n1+1):(n-h+1));

%historical mean
y_hm=zeros(n2-h+1,1);
for i=1:length(y_hm);
    y_hm(i)=mean(y(1:(n1-h+i)));  %expanding window
end;
MSE_hm=mean((y_hm-y_true).^2)
MSE_0=mean((0-y_true).^2);


%direct method
%y_direct=zeros(n2-h+1,1);
%for i=1:length(y_direct);
 %   res=ols(r(2:(n1-h+i)),[ones(n1-h+i-1,1),x(1:(n1-h+i-1))]);
  %  a_lh=res.beta(1);
  %  b_lh = res.beta(2);
   % y_direct(i)=a_lh+b_lh*x(n1-1+i);  
%end;
%MSE=mean((y_direct-y_true).^2)

%Moving Average Method 
y_MA=zeros(n2-h+1,1);
type = 'linear';
y_MA = movavg(y((n1+1):(n-h+1)),type,100);
MSE=mean((y_MA-y_true).^2)

%Moving Average Method - Lagging
y_MA_lag=zeros(n2-h+1,1);
type = 'linear';
y_MA_lag = movavg(y((n1+1):(n-h+1)),type,300);
MSE_lag=mean((y_MA_lag-y_true).^2)

%Moving Average Method - Leading 
y_MA_lead=zeros(n2-h+1,1);
type = 'linear';
y_MA_lead = movavg(y((n1+1):(n-h+1)),type,50);
MSE_lead=mean((y_MA_lead-y_true).^2)

%These are for plotting purpose
yplot_true=[y(1:n1);y_true];  
yplot_hm=[y(1:n1);y_hm];
yplot_MA=[y(1:n1);y_MA];
yplot_MA_lag=[y(1:n1);y_MA_lag];
yplot_MA_lead=[y(1:n1);y_MA_lead];

plot(date,yplot_true)
hold on 
plot(date,yplot_hm)
plot(date,yplot_MA)
plot(date,yplot_MA_lag)
plot(date,yplot_MA_lead)
hold off
legend('Actual', 'Historical Mean','Moving Average', 'Moving Average Lagging', 'Moving Average Leading');

%R-square
R2=1-MSE/MSE_hm
R2_lag = 1 - MSE_lag/MSE_hm
R2_lead = 1- MSE_lead/MSE_hm