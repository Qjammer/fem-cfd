classdef sparsetensor
	
	properties
		shape
		data
		pos
		dend
		typeof = 'sparsetensor'
	end
	
	methods
		function obj = sparsetensor(shape,allspace)
			
			if (~isvector(shape))
				error("Please, i'm not that advanced yet (Shape specified is not a vector");
			end
			
			if (size(shape,2)>1)
				shape = shape';
			end
			
			order = length(shape);
			
			if nargin < 2
				allspace = 10;
			end
			
			obj.data = zeros(allspace,1);
			obj.pos = zeros(allspace,order);
			
			obj.dend = 0;
			
			obj.shape = shape;
			
		end
		
		function o = order(self)
			o = size(self.shape,1);
		end
		
		function self = append(self,value,index)
			self.dend = self.dend+1;
			
			self.data(self.dend) = value;
			self.pos(self.dend,:) = index;
		end
		
		function [val,fvect] = at(self,index)
			n = ndims(index);
			
			if (n > self.order)
				error('Interdimensional travel forbidden (Index order greater than tensor order)');
			end
			
			if (n < self.order)
				index = [index,ones(1,self.order)];
			end
			
 			fvect = 1:size(self.pos,1);
			
			for i=1:self.order
% 				fvect(fvect) = self.pos(fvect,i)==index(i);
				fvect = fvect(self.pos(fvect,i)==index(i));
				if (isempty(self.pos))
					val = 0;
					return
				end
			end
			
			val = self.data(fvect);
		end
		
		function self = writeat(self,value,index)
			[~,f] = self.at(index);
			
			if (~any(f))
				self.append(value,index);
			else
				self.data(f) = value;
			end
		end
		
		function in = end(self,k,n)
			if (n ~= self.order)
				error('Interdimensional travel forbidden (Index order greater than tensor order)')
			end
			
			in = self.shape(k);
		end
		
		function varargout = subsref(self,in)
			
			if(in(1).type=='.') % Assures that .* calls a function, not tries to access like A(*)
				if size(in,2)<2
					varargout = self.(in(1).subs); % Accessing property
				else
					varargout = cell(1,nargout(['sparsetensor>sparsetensor.' in(1).subs]));
					[varargout{:}] = self.(in(1).subs)(in(2:end).subs{:}); % Accessing method
				end
				
				return
			end
			
			ndim = length(in.subs);
			ivector = ones(ndim,1);
			ilim = ones(ndim,1);
			
			for dim=1:ndim
				ilim(dim) = length(in.subs{dim});
			end
			
			finished = false;
			
			allmemory = ceil(size(self.data,1)*(prod(ilim)/prod(self.shape)));
			varargout = sparsetensor(ilim,allmemory);
			
			while ~finished
				
								
				
				% Updating index vector
				finished = true;
				for k=1:ndim
					ivector(k) = ivector(k) + 1;
					if (ivector(k) <= ilim(k))
						finished = false;
						break;
					end
					ivector(k) = 1;
				end
			end
		end
	end
end







