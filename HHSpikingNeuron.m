% Adapted from: https://www.math.mcgill.ca/gantumur/docs/reps/RyanSicilianoHH.pdf

classdef HHSpikingNeuron
    properties
        % Constants %%%%%%%%%%%%%%%%%
        % Reversal Potential
        SODIUM_EQUILIBRIUM{}
        POTASSIUM_EQUILIBRIUM{}
        LEAKAGE_EQUILIBRIUM{}
        
        % Capacitance
        MEMBRANE_CAPACITANCE{}
        
        % Current
        APPLIED_CURRENT{}
        
        % Conductance
        LEAKAGE_CONDUCTANCE{}
        
        % Variables %%%%%%%%%%%%%%%%%
        % Probabilities
        PotassiumGateActivation {} % n, potassium gate activation probability
        SodiumGateActivation    {} % m, sodium gate activation probability
        SodiumGateInactivation  {} % h, sodium gate inactivation probability
        
        % Conductance
        SodiumConductance{}
        PotassiumConductance{}
        
        % Current
        Amperage{}
        
        % Downstream Neurons
        Terminals{}
    end
    methods
        % Initialization
        function neuron = HHSpikingNeuron(neurons)
            % Constants
            neuron.SODIUM_EQUILIBRIUM       = 55.17;    % mv Na reversal potential
            neuron.POTASSIUM_EQUILIBRIUM    = -72.14;   % mv K reversal potential
            neuron.LEAKAGE_EQUILIBRIUM      = -49.42;   % mv Leakage reversal potential
            neuron.APPLIED_CURRENT          = 0.1;      % Applied Current
            neuron.MEMBRANE_CAPACITANCE     = 0.01;     % Membrane Capacitance
            neuron.LEAKAGE_CONDUCTANCE      = 0.003;    % mS/cm^2 Leakage conductance
            
            % Variables
            neuron.SodiumConductance        = 1.2;      % mS/cm^2 Na conductance
            neuron.PotassiumConductance     = 0.36;     % mS/cm^2 K conductance
            if exist('neurons', 'var')
                neuron.Terminals = neurons;
            end
        end
        
        % Note: α is binding constant and β is unbinding constant (I think)
        % Rate Constants (αn, βn, αm, βm, αh, βh); vary with voltage, not time
        function an=alphaRateConstantPotassiumActivation(~, voltage)   % αn, potassium gate activation rate
            an=0.01*(voltage+50)/(1-exp(-(voltage+50)/10));
        end
        function bn=betaRateConstantPotassiumActivation(~, voltage)    % βn, potassium gate activation rate
            bn=0.125*exp(-(voltage+60)/80);
        end
        function am=alphaRateConstantSodiumActivation(~, voltage)      % αm, sodium gate activation rate
            am=0.1*(voltage+35)/(1-exp(-(voltage+35)/10));
        end
        function bm=betaRateConstantSodiumActivation(~, voltage)       % βm, sodium gate activation rate
            bm=4.0*exp(-0.0556*(voltage+60));
        end
        function ah=alphaRateConstantSodiumInactivation(~, voltage)    % αh, sodium gate inactivation rate
            ah=0.07*exp(-0.05*(voltage+60));
        end
        function bh=betaRateConstantSodiumInactivation(~, voltage)     % βh, sodium gate inactivation rate
            bh=1/(1+exp(-(0.1)*(voltage+30)));
        end
        
        % Gate Activation/Inactivation Probabilities (n, m, h)
        function n=potassiumGateActivation(neuron, voltage)     % n, potassium gate activation probability
            n=neuron.alphaRateConstantPotassiumActivation(voltage)/(neuron.alphaRateConstantPotassiumActivation(voltage)+neuron.betaRateConstantPotassiumActivation(voltage));
        end
        function m=sodiumGateActivation(neuron, voltage)        % m, sodium gate activation probability
            m=neuron.alphaRateConstantSodiumActivation(voltage)/(neuron.alphaRateConstantSodiumActivation(voltage)+neuron.betaRateConstantSodiumActivation(voltage));
        end
        function h=sodiumGateInactivation(neuron, voltage)      % h, sodium gate inactivation probability
            h=neuron.alphaRateConstantSodiumInactivation(voltage)/(neuron.alphaRateConstantSodiumInactivation(voltage)+neuron.betaRateConstantSodiumInactivation(voltage));
        end
        
        % Hodgkin-Huxley Neural Spike Generator
        function dydt=hodgkinHuxley(neuron, ~, input)
            voltage = input(1);                         % V
            neuron.PotassiumGateActivation = input(2);  % n, potassium gate activation probability
            neuron.SodiumGateActivation = input(3);     % m, sodium gate activation probability
            neuron.SodiumGateInactivation = input(4);   % h, sodium gate inactivation probability
            
            % Conductance (note, leakage conductance (GL) is constant)
            neuron.SodiumConductance      = neuron.SodiumConductance    *neuron.SodiumGateActivation^3 *neuron.SodiumGateInactivation;  % GNa
            neuron.PotassiumConductance   = neuron.PotassiumConductance *neuron.PotassiumGateActivation^4;                              % GK

            % Current
            neuron.Amperage = [
                (neuron.SodiumConductance    *(voltage-neuron.SODIUM_EQUILIBRIUM));     % INa   Sodium current
                (neuron.PotassiumConductance *(voltage-neuron.POTASSIUM_EQUILIBRIUM));  % IK    Potassium current
                (neuron.LEAKAGE_CONDUCTANCE  *(voltage-neuron.LEAKAGE_EQUILIBRIUM));    % IL    Leakage current
            ];
            
            % Hodgkin-Huxley Model Differential
            dydt = [
                ((1/neuron.MEMBRANE_CAPACITANCE)                         *(neuron.APPLIED_CURRENT                -sum(neuron.Amperage)));                                                                   % dv/dt
                (neuron.alphaRateConstantPotassiumActivation(voltage)    *(1-neuron.PotassiumGateActivation)     -neuron.betaRateConstantPotassiumActivation(voltage)    *neuron.PotassiumGateActivation);  % dn/dt
                (neuron.alphaRateConstantSodiumActivation(voltage)       *(1-neuron.SodiumGateActivation)        -neuron.betaRateConstantSodiumActivation(voltage)       *neuron.SodiumGateActivation);     % dm/dt
                (neuron.alphaRateConstantSodiumInactivation(voltage)     *(1-neuron.SodiumGateInactivation)      -neuron.betaRateConstantSodiumInactivation(voltage)     *neuron.SodiumGateInactivation);   % dh/dt
            ];
        end
        
        % Solver
        function [neuralTime, neuralOutput]=spike(neuron, timespan, neuralInput)
%             disp(neuralInput)
            input = neuralInput(:,1);   % Temporary; neurons can have many
                                        % inputs. Right now we're just
                                        % using the first input.
            
            % Will be using a loop kinda like this probably
%             neuralIndex = 0;
%             while neuralIndex < length(neuralInput)
%                 neuralIndex = neuralIndex+1;
%                 
%             end
            
            voltage = input(1);                         % V
            neuron.PotassiumGateActivation = input(2);  % n, potassium gate activation probability
            neuron.SodiumGateActivation = input(3);     % m, sodium gate activation probability
            neuron.SodiumGateInactivation = input(4);   % h, sodium gate inactivation probability
            
            % If we wanted to have this function generate the activation
            % variables the following lines should be uncommented. However,
            % I'm pretty sure that breaks things at present time.
            
%             neuron.PotassiumGateActivation  = neuron.potassiumGateActivation(voltage);  % n, potassium gate activation probability
%             neuron.SodiumGateActivation     = neuron.sodiumGateActivation(voltage);     % m, sodium gate activation probability
%             neuron.SodiumGateInactivation   = neuron.sodiumGateInactivation(voltage);   % h, sodium gate inactivation probability
%             
            [time,output]=ode45(@neuron.hodgkinHuxley, timespan, [
                voltage;                            % V
                neuron.PotassiumGateActivation;     % n, potassium gate activation probability
                neuron.SodiumGateActivation;        % m, sodium gate activation probability
                neuron.SodiumGateInactivation;      % h, sodium gate inactivation probability
            ]);
            
            output = output';
            terminalIndex = 1;
            
            terminals = length(neuron.Terminals);
            
            neuralOutput = {};
            neuralTime = {};
            
            neuralOutput{terminalIndex} = output;
            neuralTime{terminalIndex} = time;
            
            % NOTE: Neurons probably should have an ID associated with them
            % so we can figure out what's what more easily. Right now it's
            % a bit of a mess.
            
            
            % TODO: Need some way of returning neurons that belong to
            % neurons; right now if a neuron object contains more neurons
            % it doesn't do anything with those additional neurons (I
            % think).
            while terminalIndex < terminals
                terminalIndex = terminalIndex+1;
                outputIndex = 0;
                while outputIndex < length(output) % TODO: need some way of integrating these inputs instead of what we're doing here
                    outputIndex = outputIndex+1;
                    [neuralTime(terminalIndex), neuralOutput(terminalIndex)] = neuron.Terminals(terminalIndex).spike(time, output(:,outputIndex));
                end
            end
        end
    end
end