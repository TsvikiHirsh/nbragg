%% Make base odf
setMTEXpref('EulerAngleConvention','Matthies');
CS = crystalSymmetry("Fm-3m",[3.3,3.3,3.3]);
ori = [orientation.cube(CS),...
    orientation.cubeND22(CS),...
    orientation.cubeND45(CS),...
    orientation.cubeRD(CS),...
    orientation.brass(CS),...
    orientation.goss(CS),...
    orientation.copper(CS)]

% Export for ncrystal
f = fopen("simple_components.csv",'w')
fprintf(f,"alpha,beta,gamma,volume,fwhm,xh,xk,xl,yh,yk,yl,zh,zk,zl\n")
% fprintf("\n========Exact==========\n")
for i=1:length(ori)
    xhkl=inv(ori(i))*xvector;
    yhkl=inv(ori(i))*yvector;
    zhkl=inv(ori(i))*zvector;
    fprintf(f,"%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n",ori(i).alpha/degree,ori(i).beta/degree,ori(i).gamma/degree,1,10,xhkl.h,xhkl.k,xhkl.l,yhkl.h,yhkl.k,yhkl.l,zhkl.h,zhkl.k,zhkl.l);
end
fclose(f);
