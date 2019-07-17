function [x_hat,B_hat]=estimateLocation(Y,lambda,Delta,W,c,useSubArray,dbar)
    % Near-Field Joint Localization and Synchronization Beyond 5G
    % (c) Henk Wymeersch, 2019
    % Usage: this code computes an estimate of the user location and clock
    % bias
    % [x_hat,B_hat]=estimateLocation(Y,lambda,Delta,W,c,useSubArray,xUE)
    % inputs:
    % Y = observation of size (K+1,N+1), for K+1 subcarriers, N+1 antennas
    % lambda = wavelength [m]
    % Delta = inter-antenna spacing [m]
    % W = bandwidth [GHz]
    % c = speed of light m/ns
    % useSubarray = 0 / 1, if we want to use the sub-array approach
    % dbar = a reference distance
        K=size(Y,1)-1;
        N=size(Y,2)-1;
        maxRange=(K+1)/W*c;        
        rangeRes=c/W;
        spatialFreqRes=lambda/(Delta*(N+1)); 
        spatialFreqRange=lambda/Delta;         
        % find N to make sure we are in far-field
        %                            in narrowband
        Ntilde=min(N+1,min(ceil(sqrt(max(dbar,1)*lambda/(2*Delta^2))),ceil(10*c/(W*Delta))));                             
        if (~useSubArray)        
            Ntilde=N;        
        end                
        if (Ntilde<8) % this means large Delta, so we only rely on measurements across frequency
            RangeEstimate=zeros(1,N+1);
            for n=-N/2:N/2
                ni=n+N/2+1;
                Kfft=K*10;
                Yn=[Y(:,ni); zeros(Kfft-K-1,1)];        % pad with zeros for higher resolution 
                [~,index]=max(abs(ifft(Yn)));           % find peak
                loc(ni,:)=[ni*Delta-N/2*Delta 0];
                RangeEstimate(ni)=mod((index-1)*rangeRes*(K+1)/(Kfft),maxRange);	% map peak to RangeEstimate
                if (RangeEstimate(ni)>maxRange/2)
                    RangeEstimate(ni)=RangeEstimate(ni)-maxRange;
                end
            end
            % now we have estimates of distance + bias
            % fill in code to perform estimation --> still to be written,
            % simple TDOA estimation
            x_hat=[];
            B_hat=[];
        else
             % process blocks of size Ntilde
            Ntest=floor((N+1)/Ntilde);
            Ntest=max(Ntest,1);
            if Ntest==1
                Ntilde=N+1;
            end            
            for n=1:Ntest
                starti=(n-1)*Ntilde+1; % first antenna of sub-array
                endi=n*Ntilde;         % last antenna of sub-array
                loc(n,:)=[(starti+Ntilde/2)*Delta-N/2*Delta 0]; % approximate array center of sub-array
                Nfft=10*N;             % let's use a large FFT in spatial domain
                Kfft=K;                                 
                Yn=Y(:,starti:endi);    % grab the observation at subarray
                Yn=[Yn zeros(K+1,Nfft-Ntilde); zeros(Kfft-K-1,Nfft)];   % pad with zeros to do larger FFT
                tmp=abs(ifft2(Yn));        
                mv=max(tmp(:));                
                [Ri,Ai]=find(tmp==mv);
                RangeEstimate(n)=mod((Ri-1)*rangeRes*(K+1)/(Kfft),maxRange);                
                if (RangeEstimate(n)>maxRange/2)
                    RangeEstimate(n)=RangeEstimate(n)-maxRange;
                end                
                SF=lambda/(Delta*(Nfft))*Ai;        % spatial frequency
                if (SF>spatialFreqRange-1) 
                    SF=SF-spatialFreqRange; 
                end
                if (SF<-spatialFreqRange+1) 
                    SF=SF+spatialFreqRange; 
                end
                AOAestimate(n)=-acos(SF)+pi;                
            end
            
            if (Ntest==1)
                % this is the standard model. We cannot estimate clock bias
                x_hat=[RangeEstimate*cos(AOAestimate) RangeEstimate*sin(AOAestimate)];
                B_hat=[];
            else
                % we only use AOA for position and then the ranges for
                % synchronization. This method can be further developed. 
                counter=0;
                for k1=1:Ntest-1
                    for k2=k1+1:Ntest                        
                        counter=counter+1;
                        % compute the bearing lines and their intersection                        
                        x_hattmp(counter,:)=intersectLines(loc(k1,:),loc(k2,:),AOAestimate(k1),AOAestimate(k2));                        
                    end
                end
                x_hat=mean(x_hattmp);   % the average of all the intersections of bearing lines
                for k1=1:Ntest
                   b_hat(k1)=norm(x_hat-loc(k1,:))-RangeEstimate(k1);   % sub-array estimate of the clock bias
                end                
                B_hat=mean(b_hat);      % average of the estimates of the clock biases. 
            end
            
        end
        