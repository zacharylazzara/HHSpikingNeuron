% Adapted from Appendix of: https://www.math.mcgill.ca/gantumur/docs/reps/RyanSicilianoHH.pdf
clc;
clear;

neurons = [
    HHSpikingNeuron();
    HHSpikingNeuron();
    HHSpikingNeuron();
    HHSpikingNeuron();
];

neuron = HHSpikingNeuron(neurons);

%Constants set for all Methods
timestep=0.04; % Time Step ms
timearray=0:timestep:25; %Time Array ms
tspan = [0,max(timearray)];

voltage=-80;                                            % Initial Membrane voltage
y0=[[voltage;                                           % V
    neuron.potassiumGateActivation(voltage);            % n, potassium gate activation probability
    neuron.sodiumGateActivation(voltage);               % m, sodium gate activation probability
    neuron.sodiumGateInactivation(voltage)],[voltage;   % V
    neuron.potassiumGateActivation(voltage);            % n, potassium gate activation probability
    neuron.sodiumGateActivation(voltage);               % m, sodium gate activation probability
    neuron.sodiumGateInactivation(voltage)]];           % h, sodium gate inactivation probability

%Matlab's ode45 function
wait = waitbar(0, "Calculating...");
[neuralTime, neuralOutput] = neuron.spike(tspan, y0);
close(wait);


figure;
outputs = length(neuralOutput);
index = 0;
while index < outputs
    index = index+1;
    time = neuralTime{index};
    output = neuralOutput{index}';

    actionPotential=output(:,1);
    potassiumActivationProbability=output(:,2);
    sodiumActivationProbability=output(:,3);
    sodiumInactivationProbability=output(:,4);
    
    subplot(2,outputs,index);
    plot(time,actionPotential);
    title(sprintf("Neural Action Potential: %d", index));
    xlabel('Time');
    ylabel('Im');

    subplot(2,outputs,index+outputs)
    plot(time,potassiumActivationProbability);
    title(sprintf("Activation Probabilities (n, m, h): %d", index));
    xlabel('Time');
    ylabel('Probability');

    hold on;
    plot(time,sodiumActivationProbability);
    plot(time,sodiumInactivationProbability);
    hold off;
end