function stop = optim_output_GS(x,optimValues,state)
global acc_i
% outputFcn_global()
%
% OutputFun for optimizers (fminunc, fmincon etc),  saving intermediate results 
% in global variable "outputFcn_global_GSrelative" for later access. 
%
% It is not supposed for live updates during optimization but for 
% later inspection, which is much more efficient. 
%
% Usage 
%   options = optimoptions( ... , 'OutputFcn',@outputFcn_global ); 
%   [XOpt,fval,exitflag,output] = fminunc(@fun, X0, options); 
%   outputFcn_global_GSrelative(k).x 
%
% See also the supplied example file. 
%
% Last Changes
%   Daniel Frisch, ISAS, 10.2020: created example
%   Daniel Frisch, ISAS, 11.2019: improved documentation 
% Created
%   Daniel Frisch, ISAS, 10.2019 
%

name=strcat("history_GS_updates_",acc_i, ".mat");
stop = false;
global outputFcn_global_GSrelative
switch state
  case 'init'
    outputFcn_global_GSrelative = struct();
    outputFcn_global_GSrelative.x = x;
    outputFcn_global_GSrelative.optimValues = optimValues;
    outputFcn_global_GSrelative.timerVal = tic;
    save(name,'outputFcn_global_GSrelative')
  case 'iter'
    ind = length(outputFcn_global_GSrelative)+1;
    outputFcn_global_GSrelative(ind).x = x;
    outputFcn_global_GSrelative(ind).optimValues = optimValues;
    outputFcn_global_GSrelative(ind).timerVal = toc(outputFcn_global_GSrelative(1).timerVal);
    save(name,'outputFcn_global_GSrelative')
  case 'done'
    ind = length(outputFcn_global_GSrelative)+1;
    outputFcn_global_GSrelative(ind).x = x;
    outputFcn_global_GSrelative(ind).optimValues = optimValues;
    outputFcn_global_GSrelative(ind).timerVal = toc(outputFcn_global_GSrelative(1).timerVal);
    save(name,'outputFcn_global_GSrelative')
  otherwise
    error('wrong switch')
end
end



