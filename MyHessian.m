% SumOfHessian(aRecon, method) : calculates the sum of the hessian matrix over all the pixels using one of the following methods:
% aRecon: Array (2d or 3d) over which to calculate the sum of the hessian
% method: Method to use. Options are 'dipHessian', '3point', 'fft', 'cyclic'

%***************************************************************************
%   Copyright (C) 2008-2009 by Rainer Heintzmann                          *
%   heintzmann@gmail.com                                                  *
%                                                                         *
%   This program is free software; you can redistribute it and/or modify  *
%   it under the terms of the GNU General Public License as published by  *
%   the Free Software Foundation; Version 2 of the License.               *
%                                                                         *
%   This program is distributed in the hope that it will be useful,       *
%   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
%   GNU General Public License for more details.                          *
%                                                                         *
%   You should have received a copy of the GNU General Public License     *
%   along with this program; if not, write to the                         *
%   Free Software Foundation, Inc.,                                       *
%   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
%**************************************************************************
%

function H=MyHessian(aRecon, method)        
s1=1/1.6;s2=0.1727;
switch method
    case '3point'
        if ndims(aRecon) == 2 || (size(aRecon,3) == 1)
            H{1,1}=dip_convolve1d(aRecon,s1*[1 -2 1],0,1); % second X derivative
            H{2,2}=dip_convolve1d(aRecon,s1*[1 -2 1],1,1); % second Y derivative
            H2a=dip_convolve1d(aRecon,[-1 0 1],0,1); %
            H{1,2}=dip_convolve1d(H2a,s2*[-1 0 1],1,1); % mixed derivative XY
        else
            H{1,1}=dip_convolve1d(aRecon,s1*[1 -2 1],0,1); % second X derivative
            H{2,2}=dip_convolve1d(aRecon,s1*[1 -2 1],1,1); % second Y derivative
            H{3,3}=dip_convolve1d(aRecon,s1*[1 -2 1],2,1); % second Y derivative
            H2a=dip_convolve1d(aRecon,[-1 0 1],0,1); %
            H{1,2}=dip_convolve1d(H2a,s2*[-1 0 1],1,1); % mixed derivative XY
            H2b=dip_convolve1d(aRecon,[-1 0 1],0,1); %
            H{1,3}=dip_convolve1d(H2b,s2*[-1 0 1],2,1); % mixed derivative XY
            H2c=dip_convolve1d(aRecon,[-1 0 1],1,1); %
            H{2,3}=dip_convolve1d(H2c,s2*[-1 0 1],2,1); % mixed derivative XY
        end
    case 'dipHessian'
        H=hessian(aRecon);
    case 'fft'
        error('fft Method is not implemented yet');
    case 'cyclic'
        error('cyclic Method is not implemented yet');
end
