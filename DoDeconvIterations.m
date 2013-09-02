% [myRes,msevalue,moreinfo,myoutput]=DoDeconvIterations(Update,startVec,NumIter)  : Performs a number of deconv iterations with a specific update rule
% Update : The update rule. Choices are:
%       'RL' : Ritchardson Lucy
%       'RLL' : Ritchardson Lucy, with a line search along each gradient direction
%       'RLF' : Ritchardson Lucy, with overrelation parameters chosen from a fixed table
%       others : all other choices are handed over to the minfuc routine (e.g. use 'lbfgs')

function [myRes,msevalue,moreinfo,myoutput]=DoDeconvIterations(Update,startVec,NumIter)
global ForcePos;
global ToEstimate;

if ForcePos && (isempty(ToEstimate) || ToEstimate==0)
    startVec=sqrt(abs(startVec)); % from now on the auxilary function is estimated and obj is the abs sqr of it.
end

moreinfo='';
if isempty(Update)
    Update='lbfgs';
end

switch Update
%     case 'K'
%     if useCuda
%         error('Khoros does not work with cuda. Please run again without the cuda flag');
%     else
%         myRes=kDeconvAxial({myim,psf},{'n',NumIter});
%     end
    case 'RL'
    myRes=startVec;
    eps=1e-6;
    for n=1:NumIter
        %[val,mygrad]=GenericErrorAndDeriv(myRes);
        [msevalue,mygrad]=GenericErrorAndDeriv(myRes);
        % reshape(dip_image(mygrad,'single'),[size(myim,1) size(myim,2) size(myim,3)])
        myRes=myRes .* (1-mygrad);
        myRes(myRes<eps)=eps;
        if msevalue ~= 0
            fprintf('Iteration %d, Richardson-Lucy iteration: Val %g\n',n,msevalue);
        else
            fprintf('Iteration %d, Richardson-Lucy iteration\n',n);
        end
        myoutput.trace.fval(n)=msevalue;
    end
    case 'RLF'  % fixed lines search steps in RL
        startSteps=[1 1];
        midSteps=[1 32 1 2 1 4 1 2 1];  % will be repeated
        endSteps=[1 1 1 1];  % to relax
        myRes=startVec;
        eps=1e-6;
        NumMidSteps=NumIter-length(startSteps)-length(endSteps);
        fullAlphaVec=repmat(midSteps,[1 floor(NumMidSteps/length(midSteps))]);
        fullAlphaVec(end+1:NumMidSteps)=midSteps(1:NumMidSteps-length(fullAlphaVec));
        fullAlphaVec=[startSteps fullAlphaVec endSteps];
        for n=1:NumIter
            [msevalue,mygrad]=GenericErrorAndDeriv(myRes);
            alpha=fullAlphaVec(n);
            myRes=myRes .* (1-alpha*mygrad);
            myRes(myRes<eps)=eps;
            if msevalue ~= 0
                fprintf('Iteration %d, accelerated Richardson-Lucy iteration: Val %g, alpha %g\n',n,msevalue, alpha);
            else
                fprintf('Iteration %d, accelerated Richardson-Lucy iteration\n',n);
            end
            myoutput.trace.fval(n)=msevalue;
        end
    case 'RLL'
    myRes=startVec;
    clear startVec;   % To save some memory
    c1=1e-4; % c1 - Sufficient Decrease for Armijo condition (1e-4)
    c2=0.4;  %  c2 - Curvature Decrease for Wolfe conditions (.2 for cg methods, .9 otherwise)
    LS_interp=1;   % with LS_interp=2 the iterations come to a halt after some point
    % For the Wolfe line-search, these interpolation strategies are available ('LS_interp'):
    %   - 0 : Step Size Doubling and Bisection
    %   - 1 : Cubic interpolation/extrapolation using new function and gradient values (default)
    %   - 2 : Mixed quadratic/cubic interpolation/extrapolation
    LS_multi=0;
    % When (LS_interp = 2), the default setting of (LS_multi = 0) uses cubic interpolation,
    % while if (LS_multi = 1) it uses quartic or quintic interpolation if more than one point are available
    progTol=1e-6;
    debug=0;
    doPlot=0;
    [msevalue,mygrad]=GenericErrorAndDeriv(myRes);
    %t=5e5;  % initial step size
    t=1; % 0.005;  % initial step size
    eps=1e-5;
    n=1;
    while n<NumIter
        d= -myRes.*mygrad;  % Direction in which the RL algorithm would decent
        gtd = mygrad'*d;    % Directional derivative
        [t,msevalue,mygrad,LSfunEvals] = WolfeLineSearch(myRes,t,d,msevalue,mygrad,gtd,c1,c2,LS_interp,LS_multi,25,progTol,debug,doPlot,1,@GenericErrorAndDeriv);
        myRes=myRes+t*d;  % Perform the update
        myRes(myRes<eps)=eps;
        n=n+LSfunEvals;  % Count the line search steps towards total number of iterations
        fprintf('Iteration %d, Linesearch-Richardson-Lucy iteration: Val %g, norm(grad) %g, step length %g\n',n,msevalue,norm(mygrad),t);
        if t< 1
            t=1;
        end
        myoutput.trace.fval(n)=msevalue;
    end
    clear d;
    clear mygrad;
    otherwise  % The Update parameter will be used as the method for the minFunc optimisation
    % using Polak Ribiere for update in the cg case
%options=struct('DerivativeCheck','off','Method','cg','Display','on','notify',1,'TolX',1e-39,'TolFun',10^-39,'MaxIter',NumIter); 
%options=struct('DerivativeCheck','off','Method','lbfgs','Display','on','verboseI',1,'notify',1,'optTol',1e-10,'progTol',1e-10,'TolX',1e-10,'TolFun',10^-29,'MaxIter',NumIter,'LS_init',3,'LS',3,'t0',1e5); 
%options=struct('DerivativeCheck','off','Method','newton0','Display','full','verboseI',1,'notify',1,'TolX',1e-10,'TolFun',10^-29,'MaxIter',NumIter,'LS_init',3,'LS',3,'t0',1e5); 
%options=struct('DerivativeCheck','off','Method','lbfgs','Display','full','verboseI',1,'notify',1,'TolX',1e-10,'TolFun',10^-29,'MaxIter',NumIter,'LS_init',3,'LS',3,'t0',1e5); 
%options=struct('DerivativeCheck','off','Method','lbfgs','Display','verbose','notify',1,'optTol',1e-10,'progTol',1e-10,'MaxIter',NumIter,'LS_type',1,'LS_interp',1,'LS_init',4,'t0',1e5); 
%options=struct('DerivativeCheck','off','Method','lbfgs','Display','verbose','notify',1,'optTol',1e-10,'progTol',1e-10,'MaxIter',NumIter,'LS_type',1,'LS_interp',1,'LS_init',4,'t0',1e5); 
%options=struct('DerivativeCheck','off','Method','lbfgs','Display','verbose','notify',1,'optTol',1e-10,'progTol',1e-10,'MaxIter',NumIter,'t0',1e5); 
%options=struct('DerivativeCheck','off','Method','lbfgs','Display',1,'verbose',1,'debug',1,'notify',1,'optTol',1e-10,'progTol',1e-10,'MaxIter',NumIter,'LS_type',1,'LS_interp',2,'t0',1e5); 
%options=struct('DerivativeCheck','off','Method','lbfgs','Display',1,'verbose',1,'debug',1,'notify',1,'optTol',1e-3,'progTol',1e-3,'MaxIter',NumIter,'MaxFunEvals',NumIter*2,'LS_type',1,'LS_interp',2,'LS_init',3,'t0',1e5); 
%options=struct('DerivativeCheck','off','Method','pcg','Display',1,'verbose',1,'debug',1,'notify',1,'optTol',1e-3,'progTol',1e-3,'MaxIter',NumIter,'MaxFunEvals',NumIter*2,'LS_type',1,'LS_interp',2,'LS_init',3,'t0',1e5); 
% options=struct('DerivativeCheck','off','Method','cg','Display',1,'verbose',1,'debug',1,'notify',1,'optTol',1e-8,'progTol',1e-8,'MaxIter',NumIter,'MaxFunEvals',NumIter*2,'LS_type',1,'LS_interp',1,'LS_init',3,'t0',1e5); 
  options=struct('useMex',0,'doPlot',1,'cgUpdate',1,'CORR',3,'DerivativeCheck','off','Method',Update,'Display',1,'verbose',1,'debug',1,'notify',1,'optTol',1e-19,'progTol',1e-19,'MaxIter',NumIter,'MaxFunEvals',NumIter*2,'LS_type',1,'LS_interp',1,'LS_init',3,'t0',1.0); 
    [myRes,msevalue,moreinfo,myoutput]=minFunc(@GenericErrorAndDeriv,startVec,options); % @ means: 'Function handle creation'     
end

if ForcePos && (isempty(ToEstimate) || ToEstimate==0)
    myRes=abssqr(myRes); % to obtain the all positive object estimate
end
