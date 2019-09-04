clear
load('/home/brian/Desktop/income.mat')

ny = income.ny;
deathrate = 1/200;
Nsim = 5e6;
Tsim = 100;

yind = zeros(Nsim,Tsim,'uint8');

yrand0 = rand(Nsim,1,'single');
yrand1 = rand(Nsim,Tsim,'single');
yrand2 = rand(Nsim,Tsim,'single');
diesim = rand(Nsim,Tsim,'single') <= deathrate;

cumdist = cumsum(income.ydist');
outrates = - diag(income.ytrans);
transrates = income.ytrans - diag(diag(income.ytrans));
transrates = transrates ./ outrates;
cumtransrates = cumsum(transrates,2);

[~,yind(:,1)] = max(yrand0<=cumdist,[],2);
for i = 1:Tsim-1
    stay = (yrand1(:,i) < 1 - outrates(yind(:,i))) & (diesim(:,i)==0);
    yind(stay,i+1) = yind(stay,i);
    
    for iy = 1:ny
        idx = (yind(:,i) == iy) & ~stay & (diesim(:,i)==0);
        [~,yind(idx,i+1)] = max(yrand2(idx,i)<=cumtransrates(iy,:),[],2);
    end
    
    idx = (diesim(:,i) == 1);
    [~,yind(idx,i+1)] = max(yrand2(idx,i)<=cumdist,[],2);
    
end

ysim = income.y.vec(yind(:,Tsim-3:Tsim));
lincome = log(sum(ysim,2));
stdev_lincome = std(lincome)