function [e] = obj2(alpha,beta,gamma,vol,halfwidth,odf)
%component_odf computes the component odf from a list of orientation and
%volume data
    odfc = (1-sum(vol))*uniformODF(odf.CS);
    ori = orientation.byEuler(alpha,beta,gamma,'ZYZ',odf.CS);
    for i=1:length(ori)
       odfc = odfc + vol(i)*unimodalODF(ori(i),SO3vonMisesFisherKernel('halfwidth',halfwidth(i)*degree));
    end
    e=calcError(odf,odfc);
end

