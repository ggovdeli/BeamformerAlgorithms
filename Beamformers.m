t = 0:0.001:0.5;                % Time range is between 0.0-0.5s, Sampling frequency is 1kHz.
s = zeros(size(t));             % Signal array is defined
s = s(:);                       % Signal in column vector
s(201:230) = s(201:230) + 1; 
s(250:290) = s(250:290) + 1;    % Defining the pulse length
figure(1);
subplot (3,1,1);
plot(t,s);                      % Plotting the pulse
title('Pulse Signal');xlabel('Time (s)');ylabel('Amplitude (V)');
carrierFrequency = 100e6;            % Carrier frequency is 100 MHz
wavelength = physconst('LightSpeed')/carrierFrequency; % Radar signal's wavelength is defined.
% Creating an ULA
ULA = phased.ULA('NumElements',10,'ElementSpacing',wavelength/2); % An ULA is defined with 10 isotropic antennas, spacing is half the wavelength.
ULA.Element.FrequencyRange = [90e5 110e6]; % Frequency range is between 90 and 110 MHz.
inputAngle = [60; 0];  % Signal comes from 60 degress azimuth and 0 degrees elevation.
collectedPW = collectPlaneWave(ULA,s,inputAngle,carrierFrequency);  % Collecting the plane wave on the ULA array.
rs = RandStream.create('mt19937ar','Seed',2008); % Random number generator for noise.
noisePower = 0.5;   % Noise power is chosen as 1 watt.
noise = sqrt(noisePower/2)*(randn(rs,size(collectedPW))+1i*randn(rs,size(collectedPW))); % Noise is defined with random generator and noise power.
receivedSignal = collectedPW + noise; % Total returned signal corresponds to received signal + noise.
subplot(3,1,2); 
plot(t,abs(receivedSignal(:,3)));axis tight;
title('Pulse Signal at Antenna 3 with Noise');xlabel('Time (s)');ylabel('Magnitude (V)'); % Antenna 3 is chosen to plot received signal.
subplot(3,1,3);
plot(t,abs(receivedSignal(:,6)));axis tight;
title('Pulse Signal at Antenna 6 with Noise');xlabel('Time (s)');ylabel('Magnitude (V)'); % Antenna 6 is chosen to plot received signal.
% Phase Shift Beamforming
PSBeamformer = phased.PhaseShiftBeamformer('SensorArray',ULA,'OperatingFrequency', ...
    carrierFrequency,'Direction',inputAngle,'WeightsOutputPort', true);  % Phase shift beamformer is created to suppress unwanted signals with a phase shift on antennas.
[yPSB,w] = PSBeamformer(receivedSignal); % Weighting coefficients and output signal is obtained.
figure(2);
subplot (3,3,1);
plot(t,abs(yPSB)); axis tight;
title('Output of Phase Shift Beamformer'); % Plotting the output
xlabel('Time (s)');ylabel('Magnitude (V)'); % Output SNR is approx. 10 times stronger (10 element array has an array gain of 10)
subplot(3,3,2);
pattern(ULA,carrierFrequency,-180:180,0,'Weights',w,'Type','powerdb','PropagationSpeed', ... 
    physconst('LightSpeed'),'Normalize',false,'CoordinateSystem','rectangular'); 
axis([-90 90 -60 0]); %Plotting the array response with weighting (this ULA has ambiguity so -90:90 degrees are plotted)
psbfLegend = legend('Phase Shift BF Response Pattern');
set(psbfLegend,'Location','best')
% Interference issue in Phase Shift Beamforming
nSamples = length(t); 
intSource1 = 10*randn(rs,nSamples,1); % Creating two interference sources 
intSource2 = 10*randn(rs,nSamples,1); % #1 interference azimuth: 45 deg, #2 interference azimuth: 15 deg
interference = collectPlaneWave(ULA,[intSource1 intSource2],[45 15; 0 0],carrierFrequency);
noisePower = 0.00001;   % Noise power, 50dB SNR to demonstrate interference's effect on phase shift beamforming.
noise = sqrt(noisePower/2)*(randn(rs,size(collectedPW))+1i*randn(rs,size(collectedPW))); % New noise value calculation
receivedSignalInt = interference + noise;     % Interferenced signal = total interference + noise.
receivedSignal = collectedPW + receivedSignalInt;             % total received signal is defined.
yPSB = PSBeamformer(receivedSignal);    % Phase shift beamformer is applied
subplot(3,3,3);
plot(t,abs(yPSB)); axis tight;    % Beamformer is plotted 
title('Output of Phase Shift BF with interference');
xlabel('Time (s)');ylabel('Magnitude (V)');
% MVDR Beamforming
MVDRBeamformer = phased.MVDRBeamformer('SensorArray',ULA,... % MVDR beamformer is defined with 60 degree signal azimuth
    'Direction',inputAngle,'OperatingFrequency',carrierFrequency,...
    'WeightsOutputPort',true);
MVDRBeamformer.TrainingInputPort = true; % Training MVDR for better results (training enabled)
[yMVDR, wMVDR] = MVDRBeamformer(receivedSignal,receivedSignalInt); % MVDR is applied to the incoming signal with interference
subplot(3,3,4);
plot(t,abs(yMVDR)); axis tight;
title('Output of MVDR BF with interference'); % MVDR Output is plotted
xlabel('Time (s)');ylabel('Magnitude (V)');
subplot(3,3,5);
pattern(ULA,carrierFrequency,-180:180,0,'Weights',wMVDR,'Type','powerdb',... % Response pattern of MVDR Beamformer
    'PropagationSpeed',physconst('LightSpeed'),'Normalize',false,...
    'CoordinateSystem','rectangular');
axis([-90 90 -80 20]);
hold on;   % for comparison with MVDR
pattern(ULA,carrierFrequency,-180:180,0,'Weights',w,... % Response pattern of Phase Shift Beamformer
    'PropagationSpeed',physconst('LightSpeed'),'Normalize',false,...
    'Type','powerdb','CoordinateSystem','rectangular');
hold off;
compLegend = legend('MVDR BF','Phase Shift BF');
set(compLegend,'Location','best')
% Self-nulling issue in MVDR Beamforming
MVDRBeamformer_selfnullDemo = phased.MVDRBeamformer('SensorArray',ULA,... % Defining an MVDR BF to illustrate self-nulling phenomenon
    'Direction',inputAngle,'OperatingFrequency',carrierFrequency,...
    'WeightsOutputPort',true,'TrainingInputPort',false);
expedtedDirection = [62; 0]; % Creating an expected signal direction that is close to target direction, creating a small signal mismatch
MVDRBeamformer_selfnullDemo.Direction = expedtedDirection; % Defining signal direction to expected direction
[ySN, wSN] = MVDRBeamformer_selfnullDemo(receivedSignal); % Plotting the newly created MVDR BF
subplot (3,3,6);
plot(t,abs(ySN)); axis tight;
title('Output of MVDR BF with Signal Direction Mismatch'); % It can be seen that the MVDR tries to suppress the target signal@60 deg because it is treated as interference and the whole beamformer fails
xlabel('Time (s)');ylabel('Magnitude (V)');               
subplot (3,3,7);
pattern(ULA,carrierFrequency,-180:180,0,'Weights',wSN,'Type','powerdb',... % Response pattern is modeled here to show that MVDR fails to find the target azimuth
    'PropagationSpeed',physconst('LightSpeed'),'Normalize',false,...
    'CoordinateSystem','rectangular');
axis([-90 90 -40 25]);
legend('MVDR BF with SDM');
% LCMV Beamformer
LCMVBeamformer = phased.LCMVBeamformer('WeightsOutputPort',true);  % Defining LCMV Beamformer
antennaSteeringVector = phased.SteeringVector('SensorArray',ULA); % Defining an antenna steering vector to prevent self-nulling
steeringAngle = antennaSteeringVector(carrierFrequency,[62 60 64]); % Steering vectors are between 60-64 degrees
LCMVBeamformer.Constraint = steeringAngle; % LCMV is constrainted to steering vectors defined above 
LCMVBeamformer.DesiredResponse = [1; 1; 1]; % Desired responses must be 1 for all directions
[yLCMV,wLCMV] = LCMVBeamformer(receivedSignal);
subplot(3,3,8);
plot(t,abs(yLCMV)); axis tight;
title('Output of LCMV BF with Signal Direction Mismatch');
xlabel('Time (s)');ylabel('Magnitude (V)');
subplot (3,3,9);
pattern(ULA,carrierFrequency,-180:180,0,'Weights',wLCMV,'Type','powerdb',... % LCMV is compared to MVDR and plotted.
    'PropagationSpeed',physconst('LightSpeed'),'Normalize',false,...
    'CoordinateSystem','rectangular');
axis([0 90 -40 35]); 
hold on; 
pattern(ULA,carrierFrequency,-180:180,0,'Weights',wSN,...
    'PropagationSpeed',physconst('LightSpeed'),'Normalize',false,...
    'Type','powerdb','CoordinateSystem','rectangular');
hold off;
compLegend2 = legend('LCMV','MVDR');
set(compLegend2,'Location','best');