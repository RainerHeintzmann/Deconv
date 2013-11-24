% [RegMat1,RegMat2]=ParseRegularisation(mycells) : Parses the parameter strings into one or two matrices
function [RegMat1,RegMat2,RegMat3]=ParseRegularisation(mycells,toReg)
global DeconvMethod
global aResampling;
%global ComplexObj;
%ComplexObj=0;
%global IntensityData;
%IntensityData=0;
%global ForcePos;
if nargin<2 
    toReg=0; % meaning object
end

NumMaxPar=13;
RegMat1 = zeros(NumMaxPar,3);
RegMat2 = zeros(NumMaxPar,3);
RegMat3 = zeros(NumMaxPar,3);
RegMat1(NumMaxPar,1)=0;  % LeastSqr as default

switch DeconvMethod
    case 'LeastSqr'
        RegMat1(10,1)=0;
        RegMat2(10,1)=0;
        RegMat3(10,1)=0;
    case 'Poisson'
        RegMat1(10,1)=1;
        RegMat2(10,1)=1;
        RegMat3(10,1)=1;
    case 'WeightedLeastSqr'
        RegMat1(10,1)=2;
        RegMat2(10,1)=2;
        RegMat3(10,1)=2;
    otherwise
        fprintf('Unknown Update method: %s\n',DeconvMethod);
        error('For update method only LeastSqr, Poisson and WeigthedLeastSqr are allowed.');
end

if isempty(mycells)
    return
end
if ~iscell(mycells)
    error('For regularisation a cell array has to be provided with the different regularisation names and arguments.')
end
if iscell(mycells{1})  % This means the user wants to parse two or three such arrays
    if numel(mycells) > 3
        error('Use one, two or three cells for regularisation parameters');
    end
    RegMat1=ParseRegularisation(mycells{1});
    if numel(mycells)>1 && ~isempty(mycells{2})
        RegMat2=ParseRegularisation(mycells{2});
    end
    if numel(mycells)>2 && ~isempty(mycells{3})
        RegMat3=ParseRegularisation(mycells{3});
    end
    return
end

for n=1:size(mycells,1)
    switch (mycells{n,1})  % This should be the token
        case 'GS'  % Args are : Lambda
            RegMat1(1,1)=mycells{n,2};
        case 'AR'  % Args are : Lambda
            if numel(mycells{n,2}) ~= 2
                error('Using AR regularisation, please provide two values in the form [lambda eps]. An eps of one is a choice. Eps > 0');
            end
            RegMat1(2,1)=mycells{n,2}(1);
            RegMat1(2,2)=mycells{n,2}(2);
        case 'TV'  % Args are : Lambda, EpsC
            if numel(mycells{n,2}) ~= 2
                error('Using TV regularisation, please provide two values in the form [lambda epsC]. An epsC of Zero means standard TV');
            end
            RegMat1(3,1)=mycells{n,2}(1);
            RegMat1(3,2)=mycells{n,2}(2);
        case 'NegSqr'  % Args are : Lambda
            RegMat1(4,1)=mycells{n,2};
        case 'GR'  % Args are : Lambda
            RegMat1(5,1)=mycells{n,2}(1);
            RegMat1(5,2)=1e-4;
            if numel(mycells{n,2}) > 1
                RegMat1(5,2)=mycells{n,2}(2);  % optional modified Good's roughness
            end
        case 'CO'  % Conchello Regularisation. Args are : Lambda
            RegMat1(11,1)=mycells{n,2};
        case 'Reuse'
            RegMat1(6,1)=1;  % Reuse what is in aRecon
        case 'StartImg'
            if ~(isa(mycells{n,2},'dip_image') || isa(mycells{n,2},'cuda'))
                error('When submitting a starting image for object or illumination, it has to be a dip_image or cuda type');
            end
            RegMat1(6,1)=1;  % Reuse what is written into aRecon below
            if toReg==0
                global aRecon;
                aRecon=mycells{n,2};
            elseif toReg==1
                global myillu;
                myillu=mycells{n,2};
            elseif toReg==2
                global myotfs;
                myotfs=mycells{n,2};
            end
        case 'Complex'
            RegMat1(7,1)=1;  
            %ComplexObj=1;
        case 'IntensityData'
            RegMat1(8,1)=1;  
            %IntensityData=1;
        case 'ForcePos'
            RegMat1(9,1)=1;  
            %ForcePos=1;
        case 'Resample'  % reconstructed object is in a different sampling than data
            aResampling=mycells{n,2};
        case 'Bg'  % include an offset intensity into the model
            RegMat1(12,1)=mycells{n,2};
        case 'Show'  % include an offset intensity into the model
            RegMat1(13,1)=1;
        otherwise
            error('For regularisation only TV, AR, GR, CO, Complex, IntensityData, ForcePos, NegSqr, Reuse, Resample, Bg and StartImg are allowed');
    end
end

