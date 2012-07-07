% This stand-alone script employs a 3-level ESN to solve a task that involves both temporal and policy abstraction
% The timescale of the lowest level is a single trial, the second level works at a scale of 20 trials and the top 
% level requires input to be maintained for 400 trials. The network can obtain performances of around 80%
% (although it may be possible to obtain perfect performance as well, by tuning the parameters; I only spent 
% about 2 hours on this task and network)
% Remco Tukker


%============== Setting the variables and helper functions

g = @(x, b)(tanh(x + b));   %activation function
minusHalf = @(x)(x - 0.5);  %helper function

netDim = 100;
topDim = 10;
midDim = 30;
bottomDim = netDim - topDim - midDim;

time = 10000;       %length of training and of testing sequence (so total length is 2x this) (800)
delta = 0.1; % 2/timesteps;      %timeconstants are multiplied by this number because of updating steps 
samplelength = time * 2;  %total length of input/output sequence

%============== Task creation, delivering sampleinput and sampleout design
%============== matrices

nrofoutputs = 8;
inputLength = 7;

%memlength top level = 400   1 input node
%memlength middle level = 20  2 input nodes
%memlength bottom level = 0   4 input nodes

sampleout = zeros(nrofoutputs,samplelength);
sampleinput = zeros(inputLength, samplelength);
       
sampleinput(1,:) = (randi(600,1,samplelength) == 1); %800
  %smaller numbers to correct for cues given close to eachother (guesstimate)
sampleinput(2,:) = (randi(30,1,samplelength) == 1); %40
sampleinput(3,:) = (randi(30,1,samplelength) == 1); %40

%make sure that the average of these inputs is 0, to avoid unwanted
%integration
sampleinput(1,:) = sampleinput(1,:) - mean(sampleinput(1,:));
sampleinput(2,:) = sampleinput(2,:) - mean(sampleinput(2,:));
sampleinput(3,:) = sampleinput(3,:) - mean(sampleinput(3,:));

sampleinput(4,:) = (randi(2,1,samplelength) == 1) - 0.5;
sampleinput(5,:) = (randi(2,1,samplelength) == 1) - 0.5;
sampleinput(6,:) = (randi(2,1,samplelength) == 1) - 0.5;
sampleinput(7,:) = (randi(2,1,samplelength) == 1) - 0.5;


lastcue  = zeros(1,3) - 1000;
for n = 1 : samplelength
    
    if (sampleinput(1,n) > 0); lastcue(1) = n; end;
    if (sampleinput(2,n) > 0); lastcue(2) = n; end;
    if (sampleinput(3,n) > 0); lastcue(3) = n; end;
    
    if ((n - lastcue(1)) < 400 )
        if ((n - lastcue(2)) < 20 )
            outnode = 1 + (sampleinput(4,n) > 0);
        else
            outnode = 3 + (sampleinput(5,n) > 0);
        end
    else
        if ((n - lastcue(3)) < 20 )
            outnode = 5 + (sampleinput(6,n) > 0);
        else
            outnode = 7 + (sampleinput(7,n) > 0);
        end
    end
    
   sampleout(outnode,n) = 5; %rest stays zero
end
    
     
%============== Setting the input weight vectors, giving us inWM connection
%============== matrix (sparse)

inputScaling = 5;  %5
inputConnectivity = 0.2; %0.1

inWM = sparse([],[],[],netDim, inputLength, round(inputConnectivity*netDim*inputLength ) );

inWM(:,1) = [ spfun( minusHalf, sprand(topDim,1,inputConnectivity)  ); zeros(midDim, 1); zeros(bottomDim, 1) ] ;    

inWM(:,2) = [ zeros(topDim,1); spfun( minusHalf, sprand(midDim,1,inputConnectivity)  ); zeros(bottomDim, 1) ] ;    
inWM(:,3) = [ zeros(topDim,1); spfun( minusHalf, sprand(midDim,1,inputConnectivity)  ); zeros(bottomDim, 1) ] ;    
                
inWM(:,4) = [ zeros(topDim,1); zeros(midDim, 1); spfun( minusHalf, sprand(bottomDim,1,inputConnectivity)  ) ] ;
inWM(:,5) = [ zeros(topDim,1); zeros(midDim, 1); spfun( minusHalf, sprand(bottomDim,1,inputConnectivity)  ) ] ;
inWM(:,6) = [ zeros(topDim,1); zeros(midDim, 1); spfun( minusHalf, sprand(bottomDim,1,inputConnectivity)  ) ] ;
inWM(:,7) = [ zeros(topDim,1); zeros(midDim, 1); spfun( minusHalf, sprand(bottomDim,1,inputConnectivity)  ) ] ;

inWM = inWM * inputScaling;

inWM1 = [inWM(1:topDim,:) ; zeros(midDim, inputLength); zeros(bottomDim, inputLength)];
inWM2 = [zeros(topDim, inputLength); inWM(topDim+1:topDim+midDim, :); zeros(bottomDim, inputLength)];
inWM3 = [zeros(topDim, inputLength); zeros(midDim, inputLength); inWM(topDim+midDim+1:netDim, :)];

%============== Setting the timeconstants and biases of the nodes in the 
%============== network, giving us biases and tau vectors

biasScaling = 0.3; %0.3

biases = randn(netDim,1) * biasScaling;

tau = [rand(topDim,1)*0.01; rand(midDim,1)*0.2; rand(bottomDim,1)*0.2+0.8 ];
tau = tau * delta ; 

%============== Create the connection matrix of the reservoir, giving us
%============== intWM (which is sparse)

topdownscaling = 1.0; 
connectivity = 0.05;

radius = 0.95;   %0.9

succeeded = 0;
while ~succeeded

    intWMtop = spfun(minusHalf, sprand(topDim, topDim, connectivity)); % left top
    intWMbottom = spfun(minusHalf,sprand(bottomDim, bottomDim, connectivity));  % right bottom
    intWMmid = spfun(minusHalf,sprand(midDim, midDim, connectivity));  % right bottom
        
    intWMbottomtop = zeros(topDim, bottomDim);  % right top
    intWMmidtop = zeros(topDim, midDim);  % right top
    intWMbottommid = zeros(midDim, bottomDim);  % right top
        
    intWMtopbottom = spfun(minusHalf,sprand(bottomDim, topDim, connectivity)) * topdownscaling;   % left bottom 
    %intWMtopbottom = zeros(bottomDim, topDim);   % left bottom 
    intWMmidbottom = spfun(minusHalf,sprand(bottomDim, midDim, connectivity)) * topdownscaling;   % left bottom 
    intWMtopmid = spfun(minusHalf,sprand(midDim, topDim, connectivity)) * topdownscaling;   % left bottom 
        
    intWM0 = [intWMtop, intWMmidtop, intWMbottomtop; intWMtopmid, intWMmid, intWMbottommid; intWMtopbottom, intWMmidbottom, intWMbottom]; 

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

timesteps = 15; 

trackmatrix = zeros(3000, netDim);
trackcount = 0;
collection = zeros(time, netDim+inputLength);
internalState = zeros(netDim,1);
oneminustau = ones(size(tau)) - tau;

for i = 1:(time * 2)
    
    if (i == (time + 1))  %training time is over, calculate output weights and start measuring performance
        outWM = (pinv(collection(:,(topDim+midDim+1):end)) * sampleout(:,1:time)')'; 
        trainoutput = outWM * collection(:,(topDim+midDim+1):end)';
        right = 0;
        wrong = 0;
        outputlist = zeros(time, nrofoutputs);
    end
    
    if ((i > 200) && (i < 401) )  %recording time
        trackcount = trackcount + 1;
        currentinput = inWM1 * sampleinput(:,i); %first give input to top half
        internalState = oneminustau.* internalState + tau.* g( (intWM*internalState), biases) + currentinput;
        trackmatrix(trackcount, :) = internalState';
        
        trackcount = trackcount + 1;
        currentinput = inWM2*sampleinput(:,i);  %then to bottom half, to make sure the distance to output doesnt influence result
        internalState = oneminustau.* internalState + tau.* g( (intWM*internalState), biases) + currentinput;
        trackmatrix(trackcount, :) = internalState';
        
        trackcount = trackcount + 1;
        currentinput = inWM3*sampleinput(:,i);  %then to bottom half, to make sure the distance to output doesnt influence result
        internalState = oneminustau.* internalState + tau.* g( (intWM*internalState), biases) + currentinput;
        trackmatrix(trackcount, :) = internalState';
        
        for j = 4:timesteps %remove input not to influence dynamics too much and let the reservoir echo
            trackcount = trackcount + 1;
            internalState = oneminustau.* internalState + tau.* g( (intWM*internalState), biases);
            trackmatrix(trackcount, :) = internalState';
        end
    
    else   %no recording time
        
        currentinput = inWM1 * sampleinput(:,i); %first give input to top 
        internalState = oneminustau.* internalState + tau.* g( (intWM*internalState), biases) + currentinput;

        currentinput = inWM2*sampleinput(:,i);  %then to mid, to make sure the distance to output doesnt influence result
        internalState = oneminustau.* internalState + tau.* g( (intWM*internalState), biases) + currentinput;
    
        currentinput = inWM3*sampleinput(:,i);  %then to mid, to make sure the distance to output doesnt influence result
        internalState = oneminustau.* internalState + tau.* g( (intWM*internalState), biases) + currentinput;
    
        for j = 4:timesteps %remove input not to influence dynamics too much and let the reservoir echo
            internalState = oneminustau.* internalState + tau.* g( (intWM*internalState), biases);
        end
        
    end
    
    totalState = [internalState;sampleinput(:,i)];
    
    if ( i < (time + 1))  %capture state for calculating weights
        collection(i,:) = totalState';
    
    else  %measure results
    
        output = outWM * totalState((topDim+midDim+1):end );
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

perf = right / (right + wrong)
meanweights = mean(mean(abs(outWM)))
testtraindiff = (trainright / (trainright + trainwrong)  ) - perf

%============== create a MDS a plot and a activity plot 

d = 3;  % number of resulting dimensions
distances = pdist(trackmatrix(50:950,:));
lowDdata = mdscale(distances, d, 'Criterion', 'strain');
figure(1);
scatter3(lowDdata(:,1), lowDdata(:,2),lowDdata(:,3));  %this only works with d = 3

figure(2);
imagesc(trackmatrix(50:950,:));
