%% Configure preferences 
clear all; close all;
setMTEXpref('EulerAngleConvention','Matthies');
pfAnnotations = @(varargin) [];
setMTEXpref('pfAnnotations',pfAnnotations);
setMTEXpref('xaxisdirection','north')
setMTEXpref('zaxisdirection','outofplane')

%% Make base odf
pf = loadPoleFigure_beartex("RPV_gamma.xpc")
CS=pf.CS; SS=pf.SS;
% figure; plot(pf{[1,2,3]}.normalize,'upper','minmax','figsize','tiny','projection','stereo','smooth','harmonicApproximation','bandwidth',15);

odf_options={'resolution',7.5*degree,'kernel',SO3vonMisesFisherKernel('halfwidth',10*degree),'iterMin',15}
odf_ini =calcODF(pf,odf_options{:});

h_plot = {Miller(1,0,0,CS),Miller(1,1,0,CS),Miller(1,1,1,CS)};
hh=figure; plotPDF(odf_ini,h_plot,'upper','minmax','figsize','tiny','projection','earea','contourf','colorRange',[0.5,2.5])

%% Extract modes 
[ori_calc, vol_calc] = calcComponents(odf_ini,'tolerance',0.5*degree);%,'exact');

[fId,ori_calc] = calcCluster(ori_calc,'maxAngle',...
       5*degree,'method','hierarchical','silent');

nori_calc = length(ori_calc);
vol = zeros(nori_calc,1);
for i=1:nori_calc
   vol(i) = sum(vol_calc(i==fId));
end
vol_calc = vol;

[vol_calc,sortIdx] = sort(vol_calc,'descend');
ori_calc = ori_calc(sortIdx);

odfc = (1-sum(vol_calc))*uniformODF(CS);
for i=1:nori_calc
   odfc = odfc + vol_calc(i)*unimodalODF(ori_calc(i),SO3vonMisesFisherKernel('halfwidth',10*degree));
end
hh=figure; plotPDF(odfc,h_plot,'upper','minmax','figsize','tiny','projection','earea','contourf','colorRange',[0.5,2.5])

%% Get halfwidths and volumes
vol = optimvar('vol',nori_calc,'LowerBound',zeros(nori_calc,1),'UpperBound',ones(nori_calc,1));
halfwidth = optimvar('halfwidth',nori_calc,'LowerBound',5*ones(nori_calc,1),'UpperBound',20*ones(nori_calc,1));
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

x0.halfwidth = (20-5)/2*ones(nori_calc,1);
x0.vol = vol_calc;

options = optimoptions('fmincon','Display','iter','MaxIterations',30,'StepTolerance',0.0001,'MaxFunctionEvaluations',5000)
[sol1,fval,exitflag] = solve(prob,x0,'Options',options)

odfc = (1-sum(sol1.vol))*uniformODF(CS);
for i=1:nori_calc
   odfc = odfc + sol1.vol(i)*unimodalODF(ori_calc(i),SO3vonMisesFisherKernel('halfwidth',sol1.halfwidth(i)*degree));
end
hh=figure; plotPDF(odfc,h_plot,'upper','minmax','figsize','tiny','projection','earea','contourf','colorRange',[0.5,2.5])

%% Try also refining orientations 
alpha = optimvar('alpha',nori_calc);
beta = optimvar('beta',nori_calc);
gamma = optimvar('gamma',nori_calc);
vol = optimvar('vol',nori_calc,'LowerBound',zeros(nori_calc,1),'UpperBound',ones(nori_calc,1));
halfwidth = optimvar('halfwidth',nori_calc,'LowerBound',5*ones(nori_calc,1),'UpperBound',20*ones(nori_calc,1));

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

hh=figure; plotPDF(odf_ini,h_plot,'upper','minmax','figsize','tiny','projection','earea','contourf','colorRange',[0.5,2.5]);

%% Export for ncrystal
f = fopen("gamma_components.csv",'w')
fprintf(f,"alpha,beta,gamma,volume,fwhm,xh,xk,xl,yh,yk,yl,zh,zk,zl\n")
% fprintf("\n========Exact==========\n")
for i=1:nori_calc
    xhkl=inv(ori(i))*xvector;
    yhkl=inv(ori(i))*yvector;
    zhkl=inv(ori(i))*zvector;
    fprintf(f,"%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n",ori(i).alpha/degree,ori(i).beta/degree,ori(i).gamma/degree,sol2.vol(i),sol2.halfwidth(i)*2,xhkl.h,xhkl.k,xhkl.l,yhkl.h,yhkl.k,yhkl.l,zhkl.h,zhkl.k,zhkl.l);
end
fclose(f);
