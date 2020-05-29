function res=ApplyHyperbolicPosModel(vec,afkt)
global savedInputHyperB;
b2 = 1.0;

res=afkt(vec);

if ~iscell(res)
    mysqrt = sqrt(b2 + abssqr(res) / 4.0);
    savedInputHyperB = {(0.5 + res / mysqrt / 4.0)};  % no abs here! This saves the whole gradient multiplyer

res= mysqrt + res / 2.0;   % This is the original simple equation, but it is numerically very unstable for small numbers!
% %     # slightly better but not good:
%     taylor1 = b2 / (2.0 * mysqrt);
%     diff = res / 2.0 + mysqrt;  % for negative values this is a difference
% %     # print('diff: ' + str(diff)+", val"+str(val)+" taylor:"+str(taylor1))
%     Order2N = res * diff;
%     mask = abs(diff / res) < 2e-4;
%     Order2N(mask) = taylor1(mask); % replace values here
%     res = taylor1 + (b2 + Order2N) / (2.0 * mysqrt);  % this should be numerically more stable
else  % For data which exists as cell array. E.g. myillu
    for n=1:size(res,2)
        mysqrt = sqrt(b2 + abssqr(res{n}) / 4.0);
        savedInputHyperB{n} = (0.5 + res{n} / mysqrt / 4.0);  % no abs here! This saves the whole gradient multiplyer

        res{n}= mysqrt + res{n} / 2.0;   % This is the original simple equation, but it is numerically very unstable for small numbers!
        
%         taylor1 = b2 / (2.0 * mysqrt);
%         diff = res{n} / 2.0 + mysqrt;  % for negative values this is a difference
%     %     # print('diff: ' + str(diff)+", val"+str(val)+" taylor:"+str(taylor1))
%     %     # if tf.abs(diff/val) < 2e-4:   % this seems a good compromise between finite subtraction and taylor series
%         Order2N = res{n} * diff;
%         mask = abs(diff / res{n}) < 2e-4;
%         Order2N(mask) = taylor1(mask); % replace values here
%         res{n} = taylor1 + (b2 + Order2N) / (2.0 * mysqrt);  % this should be numerically more stable
    end
end

