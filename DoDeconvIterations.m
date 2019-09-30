% [myRes,msevalue,moreinfo,myoutput]=DoDeconvIterations(Update,startVec,NumIter)  : Performs a number of deconv iterations with a specific update rule
% Update : The update rule. Choices are:
%       'RL' : Ritchardson Lucy
%       'RLL' : Ritchardson Lucy, with a line search along each gradient direction
%       'RLF' : Ritchardson Lucy, with overrelation parameters chosen from a fixed table
%       'RLR' : Ritchardson Lucy, with overrelation from the Biggs paper.
%       others : all other choices are handed over to the minfuc routine (e.g. use 'lbfgs')

function [myRes,msevalue,moreinfo,myoutput]=DoDeconvIterations(Update,startVec,NumIter)
% global ForcePos;
global RegularisationParameters;  % This is a matrix with all possible regularisation lambdas (and other parameters)
global allObj;
global ConvertInputToModel;
global NormFac; % To correct for the NormFactor effect on the gradient in RL

CheckOutputStop = @(x,name,i,funEvals,f,t,gtd,g,d,optCond,varargin) MarkProgress(x,name,i+1,funEvals,f,t,gtd,g,d,optCond,varargin);
CheckSimpleOutputStop = @(x,i,f)  MarkProgress(x,[],i,[],f,[],[],[],[]);

if NumIter <= 0
    [err,grad]=GenericErrorAndDeriv(startVec);
    eps = abs(NumIter);
    fprintf('Testing Gradient direction total: %g with eps = %g\n',size(startVec,1),eps);
    MaxDimensions=size(startVec,1);
    if RegularisationParameters(19,1) > 0
        MaxDimensions=RegularisationParameters(19,1);
        fprintf('Only %d out of %d dimensions are tested\n',MaxDimensions, size(startVec,1));
    end
    mygrad=grad*0;
    for d=1:MaxDimensions
        fprintf('%d ',d);
        UnitD = startVec*0;
        UnitD(d) = 1;
        mygrad(d) = (GenericErrorAndDeriv(startVec+(eps * UnitD)) - err) / eps;
        if (1)
            fprintf('Analytical: %g Numerical: %g\n',grad(d),mygrad(d));
        end
        if mod(d,40)==0
            fprintf('\n');
        end
    end
    fprintf('\n');
    if (1)
        fprintf('Imaginary components\n');
        imagOffset=floor(size(startVec,1)/2);
        if MaxDimensions < imagOffset
            for d=imagOffset+1:imagOffset+MaxDimensions+1
                fprintf('%d ',d);
                UnitD = startVec*0;
                UnitD(d) = 1;
                mygrad(d) = (GenericErrorAndDeriv(startVec+(eps * UnitD)) - err) / eps;
                if (1)
                    fprintf('Analytical: %g Numerical: %g\n',grad(d),mygrad(d));
                end
                if mod(d,40)==0
                    fprintf('\n');
                end
            end
            fprintf('\n');
        end
    end
    gradim=ConvertInputToModel(grad);  % writes result into the otfrep images
    mygradim=ConvertInputToModel(mygrad);  % writes result into the otfrep images
    
    if iscell(gradim)
        dipshow(111,cat(4,gradim{1},mygradim{1}))
    else
        dipshow(112,cat(4,gradim,mygradim))
    end
    relerror = (mygrad(1:MaxDimensions) - grad(1:MaxDimensions)) ./ mean(abs(grad(1:MaxDimensions)));
    fprintf('Max Error :%g\n',max(abs(relerror)))   % Problems are caused by the hessian operator on finite arrays at the edges
    % fprintf('Max Center Error :%g\n',max(abs(relerror(1:end-1,1:end-1))))   % Problems are caused by the hessian operator on finite arrays at the edges
    myRes=startVec;
    msevalue=err;
    moreinfo=[]; myoutput=[];
    return;
end

% if RegularisationParameters(9,1) % && (isempty(ToEstimate) || ToEstimate==0) ForcePos
%     if (isempty(ToEstimate) || ToEstimate==0)
%         startVec=sqrt(abs(startVec));
%     else
%         global myillu;
%         convertVecToIllu(startVec);  % writes the result into myillu
%         startVec=sqrt(abs(squeezeIllu())); % reads myillu and generates the 
%         startVec=convertGradToVec(startVec);
%     end
% end

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
    eps=1e-8;
    for n=1:NumIter
        %[val,mygrad]=GenericErrorAndDeriv(myRes);
        [msevalue,mygrad]=GenericErrorAndDeriv(myRes);
        if CheckSimpleOutputStop(myRes,n,msevalue)
            break;
        end

        % reshape(dip_image(mygrad,'single'),[size(myim,1) size(myim,2) size(myim,3)])
        myRes=myRes .* (1-mygrad/NormFac);
        if (RegularisationParameters(9,1)==0) % ForcePos is not selected
            myRes(myRes<eps)=eps;
        end
        if msevalue ~= 0
            fprintf('Iteration %d, Richardson-Lucy iteration: Val %g\n',n,msevalue);
        else
            fprintf('Iteration %d, Richardson-Lucy iteration\n',n);
        end
        myoutput.trace.fval(n)=msevalue;
        if iscell(allObj)
            allObj{n}=ConvertInputToModel(myRes);
        end
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
        fullAlphaVec=[startSteps fullAlphaVec endSteps] / 1.5;
        for n=1:NumIter
            [msevalue,mygrad]=GenericErrorAndDeriv(myRes);
            if CheckSimpleOutputStop(myRes,n,msevalue)
                break;
            end
            
            alpha=fullAlphaVec(n);
            myRes=myRes .* (1-alpha*mygrad/NormFac);
            if (RegularisationParameters(9,1)==0) % ForcePos is not selected
                myRes(myRes<eps)=eps;
            end
            if msevalue ~= 0
                fprintf('Iteration %d, accelerated Richardson-Lucy iteration: Val %g, alpha %g\n',n,msevalue, alpha);
            else
                fprintf('Iteration %d, accelerated Richardson-Lucy iteration. MSE is zero!\n',n);
            end
            myoutput.trace.fval(n)=msevalue;
            if iscell(allObj)
                allObj{n}=ConvertInputToModel(myRes);
            end
        end
    case 'RLR'  % RL with overrelaxation chosen as described in Bigg's paper
        Use2ndOrder=0;  % does not seem to help for longer iteration times
        veryOldRes=0;
        h=0.0;
        mygrad=0;oldgrad=0;
        myRes=startVec;
        eps=1e-6;
        for n=1:NumIter
            oldRes=myRes;
            if (n>2)
                alpha=(mygrad.' * oldgrad)./(oldgrad.' * oldgrad);
            else
                alpha=1.0;
            end
            if (n>3 && Use2ndOrder)
                alpha=sqrt((mygrad.' * mygrad)./(oldgrad.' * oldgrad));
            end
            alpha=max(min(alpha,1),eps);  % limit to between eps and 1
            if (n>3 && Use2ndOrder)
                myY=myRes+alpha*h + (alpha^2/2)*h2;
            else
                myY=myRes+alpha*h;
            end
            oldgrad= mygrad;
            clear mygrad;
            [msevalue,mygrad]=GenericErrorAndDeriv(myY);
            if CheckSimpleOutputStop(myY,n,msevalue)
                break;
            end
            mygrad= - mygrad.*myRes./NormFac;
            myRes = myY + mygrad;
            clear myY;
            if (RegularisationParameters(9,1)==0) % ForcePos is not selected
                myRes(myRes<eps)=eps;
            end
            h=myRes-oldRes;
            if (n>2 && Use2ndOrder)
                h2=myRes-2*oldRes+veryOldRes;
            end
            
            if (Use2ndOrder)
                veryOldRes = oldRes;
            else
                clear oldRes;
            end
            if msevalue ~= 0
                fprintf('Iteration %d, accelerated Richardson-Lucy iteration: Val %g, alpha %g\n',n,msevalue, alpha);
            else
                fprintf('Iteration %d, accelerated Richardson-Lucy iteration. MSE is zero!\n',n);
            end
            myoutput.trace.fval(n)=msevalue;
            if iscell(allObj)
                allObj{n}=ConvertInputToModel(myRes);
            end
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
    n2=1;
    minstep=1e-7;
    %mygrad=mygrad/NormFac; % To correct for the effect of NormFac on the gradient
    
    while n<=NumIter
        d= -myRes.*mygrad;  % Direction in which the RL algorithm would decent
        gtd = mygrad'*d;    % Directional derivative
        [t,msevalue,mygrad,LSfunEvals] = WolfeLineSearch(myRes,t,d,msevalue,mygrad,gtd,c1,c2,LS_interp,LS_multi,25,progTol,debug,doPlot,1,@GenericErrorAndDeriv);
        if CheckSimpleOutputStop(myRes,n,msevalue)
            break;
        end
        %mygrad=mygrad/NormFac; % To correct for the effect of NormFac on the gradient
        myRes=myRes+t*d;  % Perform the update
        if (RegularisationParameters(9,1)==0) % ForcePos is not selected
            myRes(dip_image(myRes<eps))=eps;  % the cast is needed for the binary adressing to work in cuda. FIX THIS
        end
        if iscell(allObj)
            allObj{n2}=ConvertInputToModel(myRes);
        end
        n=n+LSfunEvals;  % Count the line search steps towards total number of iterations
        fprintf('Iteration %d, Linesearch-Richardson-Lucy iteration: Val %g, norm(grad) %g, step length %g\n',n,msevalue,norm(mygrad),t);
        n2=n2+1;
        myoutput.trace.fval(n2)=msevalue;
        if norm(t) < minstep
            fprintf('Steplength below lower limit %g. Stopping to iterate.\n',minstep);
            break;
        end
        if t< 1
            t=1;
        end
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
  options=struct('useMex',0,'doPlot',1,'cgUpdate',1,'CORR',3,'DerivativeCheck','off','Method',Update,'Display',1,'verbose',1,'debug',1,'notify',1,'optTol',1e-19,'progTol',1e-19,'MaxIter',NumIter,'MaxFunEvals',NumIter*2,'LS_type',1,'LS_interp',1,'LS_init',3,'t0',1.0,'outputFcn',CheckOutputStop); 
  [myRes,msevalue,moreinfo,myoutput]=minFunc(@GenericErrorAndDeriv,startVec,options); % @ means: 'Function handle creation'     
end

%if RegularisationParameters(9,1) % && (isempty(ToEstimate) || ToEstimate==0)  % ForcePos
%    if (isempty(ToEstimate) || ToEstimate==0)
%        myRes=abssqr(myRes); % to obtain the all positive object estimate
%    else
%        global myillu;
%        convertVecToIllu(myRes);  % writes the result into myillu
%        myRes=abssqr(squeezeIllu()); % to obtain the all positive illumination estimate
%        myRes=convertGradToVec(myRes);
%    end
%end

function stopped=MarkProgress(x,name,i,funEvals,f,t,gtd,g,d,optCond,varargin)
global ProgressLoss;
global ProgressLossFig;
global DoStop;
global NormFac;
global ToEstimate;
global Recons;
global ConvertInputToModel;
global DeconvMask;  % since only this area is finally returned..

global RefObject;      % For checking the progress to the ground truth during iterations.  RH 2017
global RefObject_SSQ;  % For checking the progress to the ground truth during iterations.  RH 2017
global RefObject_SAbs; % For checking the progress to the ground truth during iterations.  RH 2017
global RefObject_SSQ2;  % For checking the progress to the ground truth during iterations.  RH 2017
global RefObject_SAbs2; % For checking the progress to the ground truth during iterations.  RH 2017
global RefObject_NCC; % normalized Cross correlation

drawnow(); % also needed for the DoStop user interaction in old Matlab systems
stopped=0;
if ~isempty(DoStop) && DoStop 
    stopped =1;
    return; % stop the iterations (even in minFunc)
end
if isnan(f)
    txt=sprintf('NaN appeared at iteration %d. Terminating iterations.',i);
    msgbox(txt)
    stopped =1;
    return;
end

if (size(RefObject,2)==1)
    RefObject=ConvertInputToModel(RefObject);
end
x=ConvertInputToModel(x);
if ~isempty(DeconvMask)
    x=x(DeconvMask);
    y=RefObject(DeconvMask);
else
    x=x(:);
    y=RefObject(:);
end

% Protocol the comparison to the ground truth if wanted
if ~isempty(RefObject) && ~isempty(x)
    tmp1=abs(y - x);
    RefObject_SSQ(i)=sqrt(mean(mean( tmp1.^2)));
    RefObject_SAbs(i)=mean(mean( tmp1));
    tmp2=abs(y.^2 - x.^2);
    RefObject_SSQ2(i)=sqrt(mean(mean( tmp2.^2)));
    RefObject_SAbs2(i)=mean(mean( tmp2));

    RefObject_NCC(i)=mean((y-mean(y)).*(x-mean(x)))/sqrt(var(y).*var(x));
    clear tmp1;clear tmp2;
end

if ~isempty(ProgressLossFig) && ~isempty(ProgressLossFig)
        if isempty(NormFac)
            myFactor=1.0 / prod(size(Recons));
        else
            myFactor=1.0/NormFac/ prod(size(Recons));
        end
        ProgressLoss(length(ProgressLoss)+1) = myFactor * f; % current function value
        if ~isempty(ProgressLossFig)
            figure(ProgressLossFig)
            switch ToEstimate
                case 0
                    LineStyle = 'bo-';
                case 1
                    LineStyle = 'mo-';
                case 2
                    LineStyle = 'ro-';
            end
            lP = length(ProgressLoss);
            if lP == 1
                plot(lP,ProgressLoss(end),LineStyle); % /ProgressLoss(1)
            else
                plot(lP-1:lP,ProgressLoss(end-1:end),LineStyle); % /ProgressLoss(1)
            end
            hold on;
            title('Deconvolution Progress'); xlabel('Iteration no.'); ylabel('Loss Value');
            drawnow();
        end
end

