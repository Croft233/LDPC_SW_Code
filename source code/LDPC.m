
classdef LDPC < code
    properties
    end
    
    properties (Hidden)
        Gamma   %matrix for encoding
        Phi     %matrix for encoding
        Phiinv  %inverse of Phi
        E       %part of H
        
        decoderObj  %LDPC decoder object from Communications Toolbox
        
        %Tinv %TEMPORARY DELETE
        %E    %TEMP DEL
        
        %encoding matrices
        A       %Encoding A
        B       %Encoding B
        F       %Encoding F
        pi      %
        piinv  %
        Tinv    %temporary for syndrome encoding
    end
    
    methods
    end
    
    methods (Static)
        classVerify;
        H = construct802p11n(n,k);
        H = construct802p16e(n,k);
        [Lout, Lq] = mxbpstep(H0,sign_Lq,abs_Lq);
        Z = sparse2qc(H,Z);
        
        function thispath = fullPath
            %LDPC.FULLPATH  Absolute path to the @LDPC class
            
            thispath = mfilename('fullpath');
            ind = strfind(thispath,mfilename);
            thispath = [thispath(1:ind(end)-1)];
        end
        
        
        function H = triangular(H)
            %this triangular function is more of an experiment, it doesn't
            %really do much
            [r,n] = size(H);
            for rr = r:-1:2
                cc = n-r+rr;
                
                %find a row with a one in column
                t = find(H(:,cc) == 1,1);
                assert(~isempty(t));
                H([t rr],:) = H([rr t],:);
                
                for jj = rr-1:-1:1
                    if H(jj,cc) == 1
                        H(jj,:) = mod( H(jj,:) + H(rr,:) , 2);
                    end
                end
            end
        end
        
        function [H] = matrix2alt(Hin,g)
            N = length(Hin);
            
            %sort columns
            firstOne = zeros(N,1);
            for ii = 1:N
                firstOne(ii) = find(Hin(:,ii) == 1,1,'first');
            end
            s = sort(firstOne);
            L = eye(3);L = L(:,t)
            
            [H,retry] = LDPC.local_matrix2alt(Hin,g);
            ghat = LDPC.findTheGap(H);
            if ghat ~= g && retry
                keyboard
                [H,retry] = LDPC.local_matrix2alt(H,g);
            end
            ghat = LDPC.findTheGap(H);
            assert(g == ghat,'failed to verify the desired gap');
        end
        
        function [H,retry] = local_matrix2alt(Hin,g)
            %attempt to put H into ALT form with gap g using row operations
            %can fail with retry = true, meaning the algorithm should be
            %run again on the output matrix.
            H = Hin;
            n = length(H);
            retry = false;
            
            for r = (n-g):-1:1
                c = r + g;
                
                fprintf('pivot on r=%d, c=%d\n',r,c);
                %H
                
                if H(r,c) ~= 1
                    z   = find(H(:,c) == 1);
                    upper = find(z < r,1);
                    if ~isempty(upper)
                        H([r z(upper)],:) = H([z(upper) r],:);
                    else
                        retry = true;
                        lower = find(z > r,1);
                        if isempty(lower)
                            error('All-zeros column detected')
                        else
                            H([r z(lower)],:) = H([z(lower) r],:);
                        end
                    end
                end
                %H
                assert(H(r,c) == 1)
                
                onesAbove = find(H(1:r-1,c) == 1);
                for rj = 1:length(onesAbove)
                    H(onesAbove(rj),:) = mod( H(onesAbove(rj),:) + H(r,:) , 2);
                end
                %H
                
            end
        end
        
        %use code.findTheGap instead
        %         function g = findGap(Hin)
        %             %find the "gap" of a partially lower triangular matrix
        %             %works even if the diagonal has zeros
        %             %
        %             %See also findTheGap
        %
        %             H0 = Hin;
        %             [r,n] = size(H0);
        %             k     = n-r;
        %             %u is row of last "1" in column 1
        %             v = find(H0(:,1) ~= 0,1);
        %             for ii = 2:r
        %                 %s is the contribution of column ii
        %                 s = find(H0(:,ii)  ~= 0,1) + ii - 1;
        %                 v = min(v,s);
        %             end
        %
        %             %gap g = r - v
        %             g = r - v;
        %
        %         end
        
        function qc2matlab(qc,Z)
            %QC2MATLAB   Converts quasicyclic matrix to Matlab source code
            %
            %  qc2matlab(qc,z) Converts the QC matrix to matlab-formatted
            %  code and places in on the clipboard (Mac OS only)
            %
            %  example:
            %
            %  C  = LDPC('wimax_576_0.66B.alist');
            %  qc = C.sparse2qc(24);  %circulant size must be provided
            %  LDPC.qc2matlab(qc,24); %source code on Mac OS clipboard
            %
            
            [numberOfBlocksR,numberOfBlocksN] = size(qc);
            filename = 'qc2matlab.txt';
            FID = fopen(filename,'w');
            fprintf(FID,'case %d\n',(numberOfBlocksN - numberOfBlocksR)*Z);
            fprintf(FID,'proto = [...\n');
            for rr = 1:numberOfBlocksR
                for nn = 1:numberOfBlocksN
                    fprintf(FID,'%3d ',qc(rr,nn));
                end
                if rr < numberOfBlocksR
                    fprintf(FID,';...\n');
                else
                    fprintf(FID,'];\n');
                end
            end
            fclose(FID);
            
            if ismac
                system(sprintf('pbcopy < %s',filename));
                fprintf('QC Matlab source code on your clipboard\n')
                delete(filename);
            else
                fprintf('QC Matlab source code in %s\n',filename);
            end
        end
        
        
        function H = readAlist(fname)
            % reads binary parity check matrix in "alist" format from file
            % FNAME and converts it to sparse matrix used in MATLAB routines.
            % This is an interface to matrices at http://wol.ra.phy.cam.ac.uk/mackay/codes/
            %
            % Example
            %        [H] = alist2sparse('A');   % A is the ascii file in alist format
            
            
            %   Copyright (c) 1999 by Igor Kozintsev igor@ifp.uiuc.edu
            %   $Revision: 1.1 $  $Date: 2000/03/23 $ Bug fixed by Hatim Behairy
            
            fid        = fopen(fname);
            n          = fscanf(fid,'%d',1);
            m          = fscanf(fid,'%d',1);
            maxinrow   = fscanf(fid,'%d',1);
            junk       = fscanf(fid,'%d',1); % no need
            num        = fscanf(fid,'%d',[1 n]); % number of elements in rows
            
            num2(1:n)  = maxinrow;
            junk       = fscanf(fid,'%d',[1 m]); % no need
            
            position = zeros(n,maxinrow);
            for i = 1:n
                for j = 1:num2(i)
                    position(i,j) = fscanf(fid,'%d',1);
                end
            end
            
            ii = zeros(1,sum(num));
            jj = ii;
            k = 1;
            for i=1:n
                for j=1:num(i)
                    jj(k) = i;
                    ii(k) = position(i,j);
                    ss = 1;
                    k = k+1 ;
                end
            end
            H = sparse(ii,jj,ss,m,n);
            fclose(fid);
        end
        
        
        function H = proto2sparse(proto,Z)
            H     = sparse(size(proto,1)*Z,size(proto,2)*Z);
            range = 1:Z;
            QC    = sparse(circshift(eye(Z),-1));
            
            for jj = 0:size(proto,2)-1
                
                H1 = sparse(size(proto,1)*Z,Z);
                for ii = 0:size(proto,1)-1
                    if proto(ii+1,jj+1) >= 0
                        H1(ii*Z + range,:) = QC^proto(ii+1,jj+1);
                    end
                end
                H(:,jj*Z + range) = H1;
            end
        end
        
    end
    
    
    
    
    methods
        
        function obj = LDPC(varargin)
            %OBJ = LDPC     LDPC code object
            %
            %  LDPC('filename.alist') Read LDPC code from an alist file
            %
            %  LDPC(H)  Construct from a given matrix H
            %
            %  LDPC(n,k)  Construct the corresponding code below
            %
            %  Codes from Wifi IEEE 802.11n-2009
            %
            %       n     k
            %     648   540  rate 5/6
            %     648   486  rate 3/4
            %     648   432  rate 2/3 (approx)
            %     648   324  rate 1/2
            %
            %     1296  648  rate 1/2
            %     1296  864  rate 2/3
            %     1296  972  rate 3/4
            %     1296 1080  rate 5/6
            %
            %  Codes from IEEE 802.16e-2005
            %
            %       n     k
            %     576   384  rate 2/3
            %     576   480  rate 5/6
            %     672   560  rate 5/6
            
            if nargin < 1
                H = LDPC.construct802p11n(648,324);
            elseif ischar(varargin{1})
                %read file
                H = LDPC.readAlist(varargin{1});
            elseif length(varargin{1}) == 1
                %LDPC(n,k)
                n = varargin{1};
                if nargin < 2
                    k = 0;
                else
                    k = varargin{2};
                end
                
                
                H = LDPC.construct802p11n(n,k);
                if isempty(H)
                    H = LDPC.construct802p16e(n,k);
                end
                
            elseif length(varargin{1}) > 1
                %LDPC(H)
                H = varargin{1};
            else
                error('first input must be numeric and length greater than 0');
            end
            
            obj = obj@code([],H);
            
            obj.n = size(obj.H,2);
            obj.r = size(obj.H,1);
            
            %use Matlab's Communication Toolbox decoder
            obj.decoderObj = comm.LDPCDecoder('ParityCheckMatrix',sparse(obj.H),'OutputValue','Whole codeword');
            
            %find the gap
            [g,allOnes] = LDPC.findTheGap(H);
            assert(allOnes,'Must have all ones on the diagonal')
            v = size(H,1) - g;
            
            k = obj.n - obj.r;
            
            A     = H(1:v,1:k);
            B     = H(1:v,k+1:k+g);
            T     = H(1:v,k+g+1:obj.n);
            C     = H(v+1:obj.r,1:k);
            D     = H(v+1:obj.r,k+1:k+g);
            obj.E = H(v+1:obj.r,k+g+1:obj.n);
            
            assert(all(diag(T) ~= 0),'diagonal elements must be non-zero')
            
            Tinv      = mod(inv(T),2);
            if all(unique(Tinv(:)) == [0 ; 1]) == false
                %this part may not be needed
                warning('inv(T) over reals mod 2 is non-integer')
                TT = gf(T,1);
                Tinv = inv(TT);
                Tinv = double(Tinv.x);
            end
            obj.Gamma = mod(-obj.E * Tinv * A + C,2);
            obj.Phi   = mod(-obj.E * Tinv * B + D,2);
            
            rankDeficiency  = length(obj.Phi) - rank(gf(full(obj.Phi),1));
            if rankDeficiency > 0
                error('Phi is rank deficient, you should modify the source code with a flag so that the Phi encoder cannot be used')
            else
                obj.Tinv   = mod( inv(T) , 2); %temporary until can fix encoder()
                obj.Phiinv = mod(inv(obj.Phi),2);
                obj.F      = mod( obj.Phiinv * obj.Gamma, 2);
                obj.gap    = g;
                obj.k      = obj.n - obj.r;  %dimension must be k if we can invert -- unless Phi is rank deficient
            end
        end
        
        
        function [c,u] = encoder(obj,u,s)
            %LDPC.encoder   systematic encode LDPC codeword
            %
            %   [c,u] = LDPC.encoder  randomly generates information U and
            %   encodes it to a systematic codeword C.
            %
            %   c = LDPC.encoder(u) encodes a 1-by-C.k binary vector U to a
            %   codeword C.
            %
            %   A full-rank check matrix is required.  Necessary encoding
            %   initialization is perfomred by the LDPC constructor.
            %
            
            % Phi and Gamma are special encoding matrices needed.
            if nargin < 2
                u   = randi(2,[obj.k 1]) - 1;
            else
                assert(all( size(u) == [obj.k 1] ),'input u must be k-row column vector')
            end
            
            v   = obj.r - obj.gap;
            if nargin < 3
                s       = zeros(obj.r,1);
                s1      = zeros(v,1);
                s2prime = zeros(obj.gap,1);
            else
                assert(all( size(s) == [obj.r 1] ),'input r must be r-row column vector')
                s1      = s(1:v);
                s2      = s(v+1:end);
                %Did you check this?
                %s2prime = s2 - obj.E * obj.Tinv * s1;
                s2prime = obj.Phiinv * (s2 - obj.E * obj.Tinv * s1);
                s2prime = mod(s2prime,2);
            end
            p1  = mod(- obj.F * u + s2prime,2);
            p2  = zeros(v,1);
            for ii = 1:v
                t      = obj.H(ii,:) * [u ; p1 ; p2];
                p2(ii) = mod(s1(ii) - t,2);  %assumes T(ii,ii) = 1
            end
            c = [u ; p1 ; p2];
            assert( all(mod(obj.H*c,2) == s),'Encoding error');
        end
        
        
        function [c,u] = cosetEncoder(obj,s,u)
            %LDPC.cosetEncoder  Encode to the coset of a codeword
            %
            % [c,u] = LDPC.cosetEncoder(s) Randomly generates information U
            % and encodes it to a sequence V such that H*V = S.  S is an
            % LDPC.r-by-1 vector.
            %
            % c = LDPC.cosetEncoder(u) encodes a 1-by-LDPC.k binary vector
            % U to a codeword C.
            %
            
            assert(nargin > 1,'syndrome s is required');
            assert(all(size(s) == [obj.r 1]),'Size of S must be R-by-1');
            if nargin < 3
                u   = randi(2,[obj.k 1]) - 1;
            else
                assert(all( size(u) == [obj.k 1] ),'input u must be k-row column vector')
            end
            
            [c,u] = obj.encoder(u,s);
        end
        
        
        function chat = decoder(obj,y,varChan)
            %expecting LLR inputs
            %positive LLR means x = 0
            %negative LLR means x = 1
            
            %in latticeLib, decoders inputs are generally points in space
            %construction D decoder takes mod-2 of this value.
            %still need to convert mod-2 value to LLR.
            %that is (probably) the beauty of the triangle function
            %
            LLR = (1 - 2*y) / varChan;
            
            chat = step(obj.decoderObj,LLR); %this is the Communicaitons Toolbox decoder
            %receivedBits   = step(hDec, demodSignal);
            
            
            %             hEnc = comm.LDPCEncoder;
            %             hMod = comm.PSKModulator(4, 'BitInput',true);
            %             hChan = comm.AWGNChannel(...
            %                 'NoiseMethod','Signal to noise ratio (SNR)','SNR',1);
            %             hDemod = comm.PSKDemodulator(4, 'BitOutput',true,...
            %                 'DecisionMethod','Approximate log-likelihood ratio', ...
            %                 'Variance', 1/10^(hChan.SNR/10));
            %             hDec = comm.LDPCDecoder;
            %             hError = comm.ErrorRate;
            %             for counter = 1:10
            %                 data           = logical(randi([0 1], 32400, 1));
            %                 encodedData    = step(hEnc, data);
            %                 modSignal      = step(hMod, encodedData);
            %                 receivedSignal = step(hChan, modSignal);
            %                 demodSignal    = step(hDemod, receivedSignal);
            %                 receivedBits   = step(hDec, demodSignal);
            %                 errorStats     = step(hError, data, receivedBits);
            %             end
            %             fprintf('Error rate       = %1.2f\nNumber of errors = %d\n', ...
            %                 errorStats(1), errorStats(2))
        end
        
        
        function [chat,Lout] = soriagaDecoder(obj,Lin,ldpcIter)
            %SORIAGADECODER BP LDPC Decoder suitable for turbo equalization
            %
            %  General form:
            %  [chat,Lout] = soriagaDecoder(obj,Lin,ldpcIter)
            %
            %  The soriagaDecoder is to be used in turbo equalization
            %  schemes where the LDPC decoder should maintain its internal
            %  state between turbo iterations.
            %
            %  This LDPC.soriagaDecoder is not to be used in this way,
            %  instead the source code should be used as a starting point
            %  for designing your own LDPC turbo equalization
            %  detector/decoder.  Your source code should use the function
            %  mxbpstep as LDPC.soriagaDecoder does.
            %
            %  The core source code was provided by Joseph Soriaga in
            %  December 2000.
            
            % This function is based on jldpc/mxbp.m
            
            %Make Lin sign agree with Matlab's LDPC decoder
            Lin               = -Lin;
            %            Lin               = (1 - 2*y) / varChan;
            %            Lin = y;
            H0                = obj.H;
            [~, n]            = size(H0);
            internalMessages  = zeros(1,nzmax(H0));
            [~, ~, LinSparse] = find( H0*spdiags(Lin(:),0,n,n) );
            
            % begin loop for belief propogation
            done  = 0;
            count = 1;
            
            while (count < ldpcIter + 1) && (done ~= 1)
                %  update metrics with channel prob
                if count == 1
                    internalMessages = LinSparse;
                else
                    internalMessages = internalMessages + LinSparse;
                end
                
                % do one step of bp
                [Lout, internalMessages] = LDPC.mxbpstep(H0,...
                    sign(internalMessages),abs(internalMessages));
                Lout = Lout(:);
                
                % check if codeword
                chat = (sign(Lout + Lin(:)) + 1)/2;
                if sum(mod(H0*chat(:),2))==0
                    done = 1;
                else
                    count = count + 1;
                end
            end
            
            %Make Lin sign agree with Matlab's LDPC decoder
            Lout = -Lout;
        end
        
        
        function writeAlist(obj,fname)
            %LDPC.writeAlist  Write H matrix to a file
            %
            %  C = LDPC
            %  C.writeAlist('myCode.alist');
            
            %This source code has not been debugged yet
            
            FID = fopen(fname,'w');
            assert(FID>0,'Could not open %s for writing',fname);
            fprintf('%d %d',obj.n,obj.m);
            fprintf('\n');
            columnWeights = sum(H,1);
            rowWeights = sum(H,2);
            fprintf(FID, '%d %d',max(columnWeights),max(rowWeights));
            fprintf(FID, '\n');
            fprintf(FID, '%d ',columnWeights);
            fprintf(FID, '\n');
            fprintf(FID, '%d ',rowWeights);
            for ii = 1:obj.n
                z = zeros(1,maxColumnWeights);
                t = find(obj.H(ii,:) ~= 0);
                z(1:length(t)) = t;
                fprintf(FID, '%d ',z);
                fprintf(FID, '\n');
            end
            
            for ii = 1:obj.m
                z = zeros(1,maxRowWeights);
                t = find(obj.H(:,ii) ~= 0);
                z(1:length(t)) = t;
                fprintf(FID, '%d ',z);
                fprintf(FID, '\n');
            end
            fclose(FID);
        end
        
    end
    
end
