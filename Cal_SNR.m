function SNR = Cal_SNR(G,metric,varargin)
% [ PPR PRMSR PCE ENTROPY ] = Cal_SNR( G )
% This function is used to calculate the basic SNR of the correlation plane
% The only input you need is the correlation plane (SCC,RPC)
% G is the cross-correlation plane;
% We first do the minimum correlation value subtraction to eliminate background noise effect (this will improve SCC result by a lot, but only a little for RPC)
% The math of calculating the SNR (PPR,PRMSR,PCE) can be found in Kumar's paper (1990)
%This function is written by Zhenyu Xue

CoPlane1 = G;
CoPlane = CoPlane1-min(CoPlane1(:)); %minimum subtraction to eliminate background noise effect
NX = size(CoPlane,2);
NY = size(CoPlane,1);

if isempty(varargin)
[Max(1),Ind(1)] = max(CoPlane(:)); %find the primary peak
tem = imregionalmax(CoPlane);
peakmat = CoPlane.*tem;
NPeak = 2; %number of peaks you want to include
for i = 2:NPeak
    peakmat(peakmat==Max(i-1)) = 0;
    [Max(i),Ind(i)] = max(peakmat(:)); % find the second and third peak (although the third peak is not used)
end
else % Reuse pre-computed peaks
    Max = varargin{1};
end
PEAKMAGNI = Max (1)^2;  %magnitude of the primary peak

if strcmp(metric,'PPR')
% PPR
PPR = Max (1)/Max(2);
SNR=PPR;
end

if strcmp(metric,'PRMSR')
% PRMSR
numpr = 0; % initialize number of pixels include in rms part
yy = 0; % initialize the summation of correlation value squared in the rms part
for kk = 1:NY
    for ll = 1:NX
        test = CoPlane(kk,ll);
        if test <0.5*Max(1) % pixels with correlation value less than half of primary peak will be put into the rms part
            yy = yy+test^2;
            numpr = numpr+1; %count the total number of pixels of rms part
        end
    end
end
yrms = (1/numpr*yy)^(1/2); %rms value is an average
PRMSR = PEAKMAGNI/(yrms^2);
SNR=PRMSR;
end

if strcmp(metric,'PCE')
% PCE
EE = 0; % initialize the correaltion energy (without normalized by plane size) 
for kk = 1:NY % calculate correlation energy (using the whole correaltion plane)
    for ll = 1:NX
        test2 = CoPlane(kk,ll);
        EE = EE+test2^2; % the correlation energy is the summation of the square correlation value of the whole correlation plane
    end
end
E = EE/(NY*NX); % normalized by the correlation plane size
PCE = PEAKMAGNI/E;
SNR=PCE;
end

if strcmp(metric,'ENTROPY')
% ENTROPY
for kk = 1:NY
    for ll = 1:NX
        tempCPV((ll-1)*NX+kk) = CoPlane(kk,ll); % reshape the value of the correlation plane (matrix->vector) (this can be done by using reshape)
    end
end
Nbins = 30; % # of bins to create the correlation value histogram
count = hist(tempCPV,Nbins); % bin the correlation value to genrerate the histogram
entropy = 0; % initialize entropy
for kkk = 1:Nbins
    prob(kkk) = count (kkk)/(NY*NX); % calculate the probabilty to find a certainty correlation value
    if prob (kkk) ~= 0
        entropy = entropy+prob(kkk)*log10(1/prob(kkk)); % Shannon entropy equation
    end
end
ENTROPY = entropy;
SNR=ENTROPY;
end

end
    
    