%% Configure preferences 
setMTEXpref('EulerAngleConvention','Matthies');
pfAnnotations = @(varargin) [];
setMTEXpref('pfAnnotations',pfAnnotations);
setMTEXpref('xaxisdirection','north')
setMTEXpref('zaxisdirection','outofplane')

%% Make base odf
CS = crystalSymmetry("432",[3.6,3.6,3.6]);
nori = 10
seed=1
rng(seed)
ori_ini = orientation.rand(nori,CS);
ori_ini = ori_ini.project2FundamentalRegion();
halfwidth_ub = 20;
halfwidth_lb = 5;
halfwidth_ini = randi([halfwidth_lb halfwidth_ub],nori,1)
uniform_weight_ini = 0.6;
vol = rand(10,1);
vol_ini = vol*(1-uniform_weight_ini)/sum(vol);

odfc = (1-sum(vol_ini))*uniformODF(CS);
for i=1:length(ori_ini)
   odfc = odfc + vol_ini(i)*unimodalODF(ori_ini(i),SO3vonMisesFisherKernel('halfwidth',halfwidth_ini(i)*degree));
end
odf_ini = SO3FunHarmonic(odfc.calcFourier(),CS)

h_plot = {Miller(1,0,0,CS),Miller(1,1,0,CS),Miller(1,1,1,CS)};
hh=figure; plotPDF(odf_ini,h_plot,'upper','minmax','figsize','tiny','projection','earea','contourf')

%% Extract modes 
ori_seed = equispacedSO3Grid(odf_ini.CS,odf_ini.SS,'resolution',2.5*degree)
[ori_calc, vol_calc] = calcComponents(odf_ini,'seed',ori_seed,'tolerance',1*degree);%,'exact');

nori_calc = length(ori_calc);
odfc = (1-sum(vol_calc))*uniformODF(CS);
for i=1:nori_calc
   odfc = odfc + vol_calc(i)*unimodalODF(ori_calc(i),SO3vonMisesFisherKernel('halfwidth',12*degree));
end
hh=figure; plotPDF(odfc,h_plot,'upper','minmax','figsize','tiny','projection','earea','contourf')

%% Get halfwidths and volumes
vol = optimvar('vol',nori_calc,'LowerBound',zeros(nori_calc,1),'UpperBound',ones(nori_calc,1));
halfwidth = optimvar('halfwidth',nori_calc,'LowerBound',halfwidth_lb*ones(nori_calc,1),'UpperBound',halfwidth_ub*ones(nori_calc,1));
alpha = ori_calc.alpha;
beta = ori_calc.beta;
gamma = ori_calc.gamma;

%Pass the CS to func 
func = @(halfwidth,vol)obj1(halfwidth,vol,alpha,beta,gamma,odf_ini);

%Setup the problem 
prob = optimproblem;
prob.Objective = fcn2optimexpr(func,halfwidth,vol)
prob.Constraints.cons1 = sum(vol)<=1;

show(prob)

x0.halfwidth = (halfwidth_ub-halfwidth_lb)/2*ones(nori_calc,1);
x0.vol = vol_calc;

options = optimoptions('fmincon','Display','iter','MaxIterations',30,'StepTolerance',0.0001,'MaxFunctionEvaluations',5000)
[sol1,fval,exitflag] = solve(prob,x0,'Options',options)

odfc = (1-sum(sol1.vol))*uniformODF(CS);
for i=1:nori_calc
   odfc = odfc + sol1.vol(i)*unimodalODF(ori_calc(i),SO3vonMisesFisherKernel('halfwidth',sol1.halfwidth(i)*degree));
end
hh=figure; plotPDF(odfc,h_plot,'upper','minmax','figsize','tiny','projection','earea','contourf')

%% Try also refining orientations 
alpha = optimvar('alpha',nori_calc);
beta = optimvar('beta',nori_calc);
gamma = optimvar('gamma',nori_calc);
vol = optimvar('vol',nori_calc,'LowerBound',zeros(nori_calc,1),'UpperBound',ones(nori_calc,1));
halfwidth = optimvar('halfwidth',nori_calc,'LowerBound',halfwidth_lb*ones(nori_calc,1),'UpperBound',halfwidth_ub*ones(nori_calc,1));

%Pass the CS to func 
func = @(alpha,beta,gamma,vol,halfwidth)obj2(alpha,beta,gamma,vol,halfwidth,odf_ini);
%Setup the problem 
prob = optimproblem;
prob.Objective = fcn2optimexpr(func,alpha,beta,gamma,vol,halfwidth)
prob.Constraints.cons1 = sum(vol)<=1;
show(prob)

x0.alpha = ori_calc.alpha;
x0.beta = ori_calc.beta;
x0.gamma = ori_calc.gamma;
x0.vol = sol1.vol;
x0.halfwidth = sol1.halfwidth;

options = optimoptions('fmincon','Display','iter','MaxIterations',50,'StepTolerance',0.0001,'MaxFunctionEvaluations',5000)
[sol2,fval,exitflag] = solve(prob,x0,'Options',options)

% Plot results
odfc = (1-sum(sol2.vol))*uniformODF(CS);
ori = orientation.byEuler(sol2.alpha,sol2.beta,sol2.gamma,'ZYZ',CS);
for i=1:length(ori)
   odfc = odfc + sol2.vol(i)*unimodalODF(ori(i),SO3vonMisesFisherKernel('halfwidth',sol2.halfwidth(i)*degree));
end

hh=figure; plotPDF(odfc,h_plot,'upper','minmax','figsize','tiny','projection','earea','contourf');hold on

%% Correlate extracted modes with initial modes
ind = zeros(nori,1);
for i=1:length(ori)
    ind(angle(ori(i),ori_ini) < 5*degree)=i;
end
fprintf("Missed component volume:%.2f\n",sum(vol_ini(ind==0)));

for i=1:length(ori)
    fprintf("Misorientation of Calculated Component: %0.2f\n",angle(ori(i),mean(ori_ini(ind==i)))/degree);
end

for i=1:length(ori)
    fprintf("Deviation in Volume: %0.4f\n",sol2.vol(i)-sum(vol_ini(ind==i)));
end

