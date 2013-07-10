% [RegMat1,RegMat2]=ParseRegularisation(mycells) : Parses the parameter strings into one or two matrices
function [RegMat1,RegMat2]=ParseRegularisation(mycells)
RegMat1 = zeros(5,3);
RegMat2 = zeros(5,3);
if isempty(mycells)
    return
end
if ~iscell(mycells)
    error('Dammit!! Learn how to use this program!')
end
if iscell(mycells{1})  % This means the user wants to parse two such arrays
    if numel(mycells) ~= 2
        error('Use one or two cells for regularisation parameters');
    end
    RegMat1=ParseRegularisation(mycells{1});    
    RegMat2=ParseRegularisation(mycells{2});
    return
end

for n=1:size(mycells,1)
    switch (mycells{n,1})  % This should be the token
        case 'GS'  % Args are : Lambda
            RegMat1(1,1)=mycells{n,2};
        case 'AR'  % Args are : Lambda
            RegMat1(2,1)=mycells{n,2};
        case 'TV'  % Args are : Lambda, EpsC
            if numel(mycells{n,2}) ~= 2
                error('Using TV regularisation, please provide two values in the form [lambda epsC]. An epsC of Zero means standard TV');
            end
            RegMat1(3,1)=mycells{n,2}(1);
            RegMat1(3,2)=mycells{n,2}(2);
        case 'NegSqr'  % Args are : Lambda
            RegMat1(4,1)=mycells{n,2};
        case 'GR'  % Args are : Lambda
            RegMat1(5,1)=mycells{n,2};
        otherwise
            error('For regularisation only TV, AR, GR and NegSqr are allowed');
    end
end
