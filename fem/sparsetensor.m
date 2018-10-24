classdef sparsetensor < handle
	
	properties
		shape
		data
		pos
		dend
	end
	
	methods
		function obj = tensor(shape,allspace)
			
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
		
		function val = readat(self,index)
			n = ndims(index);
			
			if (n > self.order)
				error('Interdimensional travel forbidden (Index order greater than tensor order)');
			end
			
			if (n < self.order)
				index = [index,ones(1,self.order)];
			end
			
			fvect = true(size(self.pos,1),1);
			
			for i=1:self.order
				fvect(fvect) = self.pos(fvect,i)==index(i);
				if (~any(fvect))
					break;
				end
			end
			
			val = self.data(fvect);
			
			if isempty(val)
				val = 0;
			end			
		end
	end
end





