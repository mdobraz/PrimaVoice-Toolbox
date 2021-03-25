function [xpd,csypd,ypd] = boot_pd(iterations,bootfun,data,varargin)


str = [];
if ~isempty(varargin)
	
	for i = 1:length(varargin)
		str = [str ',' sprintf('varargin{%i}',i)];
	end
end

eval(sprintf('bootdist = bootstrp(iterations,bootfun,data%s);',str));
% bootdist = bootstrp(iterations,bootfun,data);

nbins = 1000;
pd_D = fitdist(bootdist,'kernel');
xpd = min(bootdist):1/nbins:max(bootdist);
ypd = pdf(pd_D,xpd);
csypd = cumsum(ypd);
csypd = csypd ./ max(csypd);