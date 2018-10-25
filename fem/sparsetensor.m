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
				error("Please, i'm not that advanced yet (Shape specified must be a vector");
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
		
		function self = append(self,value,varargin)
			self.dend = self.dend+1;
			
			self.data(self.dend) = value;
			self.pos(self.dend,:) = cell2mat(varargin);
		end
		
		function checkIndex(self,varargin)
			
			n = length(varargin);
			
			if (n ~= self.order)
				error('Interdimensional travel forbidden (Index order greater than tensor order)');
			end
			
			for dim=1:n
				if length(varargin{dim})>1
					error("Index provided doesn't access to a single point");
				end
				
				if varargin{dim}<1 || varargin{dim}>self.shape(dim)
					error ("Index array out of bounds");
				end
			end
		end
		
		function o = out(self,index)
			
			index = self.checkIndex(index);
			
			o = false;
			
			for dim=1:length(index)
				
				if (index(dim) < 1 || index(dim) > self.shape(dim))
					o=true;
					return
				end
				
			end
			
		end
		
		function [val,fvect] = at(self,varargin)
			
			self.checkIndex(varargin{:});
			
 			fvect = 1:size(self.pos,1);
			
			for i=1:self.order
				fvect = fvect(self.pos(fvect,i)==varargin{i});
				if (isempty(fvect))
					val = 0;
					return
				end
			end
			
			val = self.data(fvect);
		end
		
		function self = writeat(self,value,varargin)
			
			% Index is a cell array
			
			[~,f] = self.at(varargin{:});
			
			if (~any(f))
				self = self.append(value,varargin{:});
			else
				self.data(f) = value;
			end
		end
		
		function self = bulkWriteAt(self,B,varargin)
			% Check dimensions
			
			ndim = length(varargin);
			
			if (ndim ~= self.order)
				error("Index array dimension mismatch");
			end
			
			if class(B)=="double" % Standard matrix here
				% Check if the matrix is provided accordingly to the index
				
				Blogor = [];
				inlogor = [];
				
				for dim=1:ndim					
					if varargin{dim} == ':'
						varargin{dim} = 1:self.shape(dim);
					end
					
					if size(B,dim)>1
						Blogor = [Blogor,dim];
					end
					
					if length(varargin{dim})>1
						inlogor = [inlogor,dim];
					end
				end
				
				% We use the dimensions where the size is greater than 1 to
				% iterate
				
				if (length(Blogor)~=length(inlogor))
					error ("Matrix being assigned isnt consistant with index");
				end
				
				Bor = length(Blogor);
				
				for dim=1:Bor					
					if size(B,Blogor(dim))~=length(varargin{inlogor(dim)})
						error ("Matrix being assigned isnt consistant with index");
					end
				end
				
				ivector = ones(Bor,1);
				ivector = num2cell(ivector);
				gvector = ones(self.order,1);
				gvector = num2cell(gvector);
				
				ilim = ones(Bor,1);
				
				for dim=1:Bor
					ilim(dim) = size(B,Blogor(dim));
					gvector{dim} = varargin{dim}(1);
				end
				
				finished = false;
				
				while ~finished
					
					% B(vector) doesnt work...
					self = self.writeat(B(ivector{:}),gvector{:});
					
					
					% Updating index vector
					finished = true;
					for k=1:Bor
						ivector{k} = ivector{k} + 1;
						if (ivector{k} <= ilim(k))
							finished = false;							
							break;
						end
						ivector{k} = 1;
					end					
					for dim=1:Bor
						gvector{inlogor(dim)} = varargin{inlogor(dim)}(ivector{Blogor(dim)});
					end
				end
				
			end
		end
		
		function sub = subtensor(self,varargin)
			
			ndim = length(varargin);
			
			if (ndim ~= self.order)
				error("Index array dimension mismatch");
			end
			
			ivector = ones(ndim,1);
			ivector = num2cell(ivector);
			gvector = ones(ndim,1);
			gvector = num2cell(gvector);
			ilim = ones(ndim,1);
			
			for dim=1:ndim
				ilim(dim) = length(varargin{dim});
				gvector{dim} = varargin{dim}(1);
			end
			
			if prod(ilim)==1 % Directly return the value, no need to create a subtensor
				sub = self.at(gvector);
				return
			end
			
			A = prod(ilim) < self.dend;
			
			
			if 1
				finished = false;

				allmemory = ceil(size(self.data,1)*(prod(ilim)/prod(self.shape)));
				sub = sparsetensor(ilim,allmemory);

				while ~finished
					
					sub = sub.subsasgn(substruct('()',ivector),self.at(gvector{:}));				
% 					sub(ivector{:}) = self.at(gvector{:});
					
					
					% Updating index vector
					finished = true;
					for k=1:ndim
						ivector{k} = ivector{k} + 1;
						if (ivector{k} <= ilim(k))
							finished = false;
							break;
						end
						ivector{k} = 1;
					end
					for k=1:ndim
						gvector{k} = varargin{k}(ivector{k});
					end
				end
			else
				
			end
						
		end
		
		function in = end(self,k,n)
			if (n ~= self.order)
				error('Interdimensional travel forbidden (End placed on a non-reachable dimension)')
			end
			
			in = self.shape(k);
		end
		
		function varargout = subsref(self,in)
			
			if(in(1).type=='.') % Assures that .* calls a function, not tries to access like A(*)
				if size(in,2)<2
					varargout{1} = self.(in(1).subs); % Accessing property
				else
					varargout = cell(1,nargout(['sparsetensor>sparsetensor.' in(1).subs{:}]));
					[varargout{:}] = self.(in(1).subs)(in(2:end).subs{:}); % Accessing method
				end
				
				return
			end
			
			if(in(1).type=="()")
				varargout{1} = self.subtensor(in(1).subs{:});
			end
						
		end
		
		function varargout = subsasgn(self,in,B)
			
			if(in(1).type=="()")
				varargout{1} = self.bulkWriteAt(B,in(1).subs{:});
				return
			end
			
		end
	end
end







