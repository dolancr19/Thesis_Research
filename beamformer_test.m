clear variables
load data42.mat
x=bandpass(data(:,1),[24500 25500],100000);
y=bandpass(data(:,2),[24500 25500],100000);
Z=[x y];
c=1500;
Fs=100000;

first=false;
start=0;
count=1;
diff=0;
for ii=1:length(Z)
    delta=mod(ii,10);
    if delta==0
        diff=peak2peak(Z(ii-9:ii,1));
        if diff>.05 && first~=true
            start=ii;
            first=true;
        end
        if diff<.01 && first==true
            stop=ii;
            first=false;
            signal=([Z(start:stop,1) Z(start:stop,2)]); 
            %Adjust ArrayNormal depending on array orientation (default is z)
            array=phased.UCA('NumElements',2,'Radius',.019,'ArrayNormal','z');
            array.Element.FrequencyRange = [20 50000];

            estimator = phased.GCCEstimator('SensorArray',array,...
                'PropagationSpeed',c,'SampleRate',Fs);%,'SensorPairSource','Property', 'SensorPair',[3,4;1,2]);
            ang(:,count) = estimator(signal);
            count=count+1;
            %release(estimator);
        
        end
    end
end

        

%d = 0.25;
%N = 5;
%mic = phased.OmnidirectionalMicrophoneElement;
%array = phased.URA([N,N],[d,d],'Element',mic);

%arrivalAng = [17;0];
%collector = phased.WidebandCollector('Sensor',array,'PropagationSpeed',c,...
%    'SampleRate',Fs,'ModulatedInput',false);
%signal = collector(y,arrivalAng);

%estimator = phased.GCCEstimator('SensorArray',array,...
%    'PropagationSpeed',c,'SampleRate',Fs);
%ang = estimator(signal)