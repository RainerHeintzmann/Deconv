function thegrad=JoinFunctions(thegrad,viewNum,fktList,numArgs)  % This is obtained using the Wirtinger Derivative of the square model
for n =1:length(fktList)
    if numArgs(n)==1
        thegrad=fktList{n}(thegrad);
    elseif numArgs(n)==2
        thegrad=fktList{n}(thegrad,viewNum);
    else
        error('unknown number of arguments')
    end
end
