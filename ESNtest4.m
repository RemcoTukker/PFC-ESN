function [perf, meanweights, testtraindiff] = ESNtest4(task, memlength, inputtype, relevantcues, inputLength, ...
    topdownscaling, netDim, connectivity, inputConnectivity, radius, inputScaling, biasScaling, timesteps, topratio)

% This function takes a large number of parameters to set up the ESN and the task. It returns the performance of the ESN, 
% the mean output weight and the performance difference between the training and test trials.
% 
% Remco Tukker
%
% parameters: 
% task: 1 policy (WCS) 2 temporal (n-back)
% memlength: 0 for no temporal abstraction, larger integers for n back task
% inputtype: 1 seperated (topological) inputs 2 uniform input
% relevantcues & inputlength: amount of relevant cues in WCS task and total number of inputs. Should match, check code
% topdownscaling: scaling factor for connections from top to bottom area in the connection matrix
% netDim: reservoir size
% connectivity: fraction of non-zero elements in the connectivity matrix
% inputConnectivity: fraction of non-zero elements in the input matrix
% radius: spectral radius, determines weight size in connectivity 
% inputScaling: input matrix is multiplied with this number
% biasScaling: how large the biases in the reservoir are
% timesteps: number of timesteps from giving input until readout
% topratio: fraction of nodes that is in the top cluster of the network

%============== Setting the variables and helper functions

g = @(x, b)(tanh(x + b));   %activation function
minusHalf = @(x)(x - 0.5);  %helper function

topDim = round(netDim*topratio);  %the dimensions of the top and bottom part of the reservoir (/4)
bottomDim = netDim - topDim;

time = 800;       %length of training and of testing sequence (so total length is 2x this) (800)
delta = 0.1; % 2/timesteps;      %timeconstants are multiplied by this number because of updating steps 
samplelength = time * 2;  %total length of input/output sequence

%============== Task creation, delivering sampleinput and sampleout design
%============== matrices

switch task
    case 1
        nrofoutputs = 2*2^relevantcues;
    case 2
        nrofoutputs = 8;
end

sampleout = zeros(nrofoutputs,samplelength);
sampleinput = zeros(inputLength, samplelength);

switch task  %create the input/output sequence
    case 1  % hierarchical
        
        for k = 1:inputLength
            sampleinput(k,:) = (randi(2,1,samplelength) == 1) - 0.5;
        end
        
        for n = 1 : samplelength
            inputnumber = 1;
            if (sampleinput(1,n) > 0)
                for k = 1:relevantcues
                    inputnumber = inputnumber + 2^(k - 1) * (sampleinput(2*k,n) + 0.5);
                end
            else
                for k = 1:relevantcues
                    inputnumber = inputnumber + 2^(k - 1) * (sampleinput(2*k+1,n) + 0.5);
                end
                inputnumber = inputnumber + 2^relevantcues;
            end
            sampleout(inputnumber,n) = 5; %rest stays zero
        end
    
    case 2  %non-hierarchical
        
        for k = 1:inputLength
            sampleinput(k,:) = (randi(2,1,samplelength) == 1) - 0.5;
        end
        
        for n = 1:samplelength
            inputnumber = 4*(sampleinput(1,n)+0.5) + 2*(sampleinput(2,n)+0.5) + (sampleinput(3,n)+0.5) + 1;  %from 1 to 8
    
            sampleout(inputnumber,n) = 5; %rest stays zero
        end
        
    otherwise
        disp('invalid task number');
end

%shift back one input in case of memlength
sampleinput(1,:)  = circshift( sampleinput(1,:) ,[0 -memlength]);
     
%============== Setting the input weight vectors, giving us inWM connection
%============== matrix (sparse)

inWM = sparse([],[],[],netDim, inputLength, round(inputConnectivity*netDim*inputLength ) );
switch inputtype
    case 1  %seperated inputs, first input to the top, rest down
        inWM(:,1) = [ spfun( minusHalf, sprand(topDim,1,inputConnectivity)  ); zeros(bottomDim, 1) ] ;    
        for k = 2:inputLength
            inWM(:,k) = [zeros(topDim, 1); spfun( minusHalf, sprand(bottomDim,1,inputConnectivity) )];
        end
    case 2  %non-seperated inputs
        inWM(:,1) = spfun(minusHalf, sprand(netDim, 1, inputConnectivity * topDim/(topDim+bottomDim) ) );
        for k = 2:inputLength
            inWM(:,k) = spfun(minusHalf, sprand(netDim, 1, inputConnectivity* bottomDim/(topDim+bottomDim) ) );
        end % seperate columns to make sure every input receives the same amount of connections    
    otherwise
        disp('invalid input type number');
end

inWM = inWM * inputScaling;

inWM1 = [inWM(1:topDim,:) ; zeros(bottomDim, inputLength)];
inWM2 = [zeros(topDim, inputLength); inWM(topDim+1:netDim, :) ];

%============== Setting the timeconstants and biases of the nodes in the 
%============== network, giving us biases and tau vectors

biases = randn(netDim,1) * biasScaling;

lowertaubound = 0.8 * exp(-0.2 * memlength); 
tau = [rand(topDim,1)*0.2+lowertaubound; rand(bottomDim,1)*0.2+0.8 ];
tau = tau * delta ; 

%============== Create the connection matrix of the reservoir, giving us
%============== intWM (which is sparse)

succeeded = 0;
while ~succeeded

    intWMtop = spfun(minusHalf, sprand(topDim, topDim, connectivity)); % left top
    intWMbottom = spfun(minusHalf,sprand(bottomDim, bottomDim, connectivity));  % right bottom
    intWMbottomtop = zeros(topDim, bottomDim);  % right top
    intWMtopbottom = spfun(minusHalf,sprand(bottomDim, topDim, connectivity)) * topdownscaling;   % left bottom 
    intWM0 = [intWMtop, intWMbottomtop; intWMtopbottom, intWMbottom]; 

    succeeded = 1;
    try
        maxval = max(abs(eigs(diag(tau)*intWM0 + (diag(ones(size(tau)))- diag(tau)))));
        intWM = intWM0 / maxval * radius;
    catch exception
        disp('Matrix creation failed, trying again');
        succeeded = 0;  %eigs failed, try a new matrix (succeeded is put to 0)
    end

end  

%==============
%============== Run the network!!
%==============

collection = zeros(time, netDim+inputLength);
internalState = zeros(netDim,1);
oneminustau = ones(size(tau)) - tau;

for i = 1:(time * 2)
    
    if (i == (time + 1))  %training time is over, calculate output weights and start measuring performance
        outWM = (pinv(collection(:,(topDim+1):netDim)) * sampleout(:,1:time)')'; 
        trainoutput = outWM * collection(:,(topDim+1):netDim)';
        right = 0;
        wrong = 0;
        outputlist = zeros(time, nrofoutputs);
    end
    
    currentinput = inWM1 * sampleinput(:,i); %first give input to top half
    internalState = oneminustau.* internalState + tau.* g( (intWM*internalState), biases) + currentinput;

    currentinput = inWM2*sampleinput(:,i);  %then to bottom half, to make sure the distance to output doesnt influence result
    internalState = oneminustau.* internalState + tau.* g( (intWM*internalState), biases) + currentinput;
    
    for j = 3:timesteps %remove input not to influence dynamics too much and let the reservoir echo
        internalState = oneminustau.* internalState + tau.* g( (intWM*internalState), biases);
    end
    
    totalState = [internalState;sampleinput(:,i)];
    
    if ( i < (time + 1))  %capture state for calculating weights
        collection(i,:) = totalState';
    
    else  %measure results
    
        output = outWM * totalState((topDim+1):netDim );
        outputlist(i - time,:) = output;
    
        target = -1;
        maxout = 0;
        maxoutnr = 0;
        for k=1:nrofoutputs
            if (sampleout(k,i) > 0.5)
               target = k; 
            end
            if (output(k) > maxout)
               maxout = output(k);
               maxoutnr = k;
            end
        end
        
        if (target == maxoutnr)
           right = right +1; 
        else
            wrong = wrong +1;
        end
    end
    
end

trainright = 0;   %measure training results
trainwrong = 0;
for i = 1:time
        target = -1;
        maxout = 0;
        maxoutnr = 0;
        for k=1:nrofoutputs
            if (sampleout(k,i) > 0.5)
               target = k; 
            end
            if (trainoutput(k,i) > maxout)
               maxout = trainoutput(k,i);
               maxoutnr = k;
            end
        end
        
        if (target == maxoutnr)
           trainright = trainright +1; 
        else
            trainwrong = trainwrong +1;
        end
end

%============== calculate output of this function

perf = right / (right + wrong);
meanweights = mean(mean(abs(outWM)));
testtraindiff = (trainright / (trainright + trainwrong)  ) - perf;

end
