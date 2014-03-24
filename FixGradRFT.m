function transformed=FixGradRFT(transformed,amask)
if ndims(transformed)<3
    transformed(:,1:end-1)=transformed(:,1:end-1)*2; % Aurelie & Rainer to make the gradient correct
else
    transformed(:,1:end-1,:)=transformed(:,1:end-1,:)*2; % Rainer
    if ndims(amask) <3 || size(amask,3) == 1
        if ndims(amask) <3
            transformed=squeeze(transformed(:,:,0));
        else
            transformed=transformed(:,:,0);
        end
    end
end
