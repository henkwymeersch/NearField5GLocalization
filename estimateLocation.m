function [x_hat,B_hat]=estimateLocation(Y,lambda,Delta,W,c,useSubArray,dbar,bias,xUE,visualize)
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
    % bias = true value of the bias, used when useSubarray=0
    % xUE = true value of the UE location, used for visualization
    % visualize = 0 / 1 to visualize some results
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
            Ntilde=N+1;        
        end                
        if (Ntilde<8) % this means large Delta, so we only rely on measurements across frequency
            RangeEstimate=zeros(1,N+1);
            for n=-N/2:N/2
                ni=n+N/2+1;
                Kfft=K+1;
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
                Nfft=N+1;                 % let's use a large FFT in spatial domain
                Kfft=K+1;                                 
                Yn=Y(:,starti:endi);    % grab the observation at subarray                
                Yn=[Yn zeros(K+1,Nfft-Ntilde); zeros(Kfft-K-1,Nfft)];   % pad with zeros to do larger FFT                
                tmp=abs(ifft2(Yn));        
                if (visualize)
                    figure(1)
                    mesh(1:Nfft,1:Kfft,tmp)
                    xlabel('spatial index')
                    ylabel('frequency index')
                    title(['# of antennas used = ' num2str(Ntilde)]);
                    pause(0.1)
                end
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
            
            if (or(Ntest==1,var(AOAestimate)==0))                
                if (useSubArray)
                    % we don't have enough information to provide the UE
                    % location
                    x_hat=[];
                    B_hat=[];                    
                else
                    % the standard approach
                    RangeEstimate=RangeEstimate+bias;
                    x_hat=[RangeEstimate*cos(AOAestimate) RangeEstimate*sin(AOAestimate)];
                    B_hat=[];                    
                end
                
            else
               % we only use AOA for position and then the ranges for
               % synchronization. This method can be further developed.   
               
               %x_hat=intersectLines(loc(1,:),loc(end,:),AOAestimate(1),AOAestimate(end)); 
                if (mod(Ntest,2)==0)                   
                    x_hat=intersectLines(loc(Ntest/2,:),loc(Ntest/2+1,:),AOAestimate(Ntest/2),AOAestimate(Ntest/2+1));                
                else
                    x_hat=intersectLines(loc(floor(Ntest/2),:),loc(floor(Ntest/2)+1,:),AOAestimate(floor(Ntest/2)),AOAestimate(floor(Ntest/2)+1));                
                end
                if (visualize) 
                    close all
                    figure(2);
                    for k1=1:Ntest
                        line([loc(k1,1) loc(k1,1)+ 4*dbar*cos(AOAestimate(k1)) ], [0 4*dbar*sin(AOAestimate(k1)) ])
                        hold on                    
                    end
                    plot(xUE(1),xUE(2),'r*')
                    plot(x_hat(1),x_hat(2),'b*')
                    hold off
                    pause(1)
                end                
                for k1=1:Ntest
                   b_hat(k1)=norm(x_hat-loc(k1,:))-RangeEstimate(k1);   % sub-array estimate of the clock bias
                end                
                B_hat=mean(b_hat);      % average of the estimates of the clock biases. 
            end
        end       