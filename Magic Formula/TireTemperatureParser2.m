function [SA_out, FY_out, FZ_out, P_out, TSTC_out, ET_out, IA_out, NFY_out, index] = TireTemperatureParser2(ET_input, SA_input, file_name)

%SA_tol = 0.5;
SA_tol = 0.5;

load(file_name);

index = [];

%isolate data based on given parameters
for i = 1:numel(ET_input)
    for j = 1:numel(SA_input)
        index = [index; find(ET < ET_input(2) & ET > ET_input(1)... 
        & SA > SA_input(j) - SA_tol & SA < SA_input(j) + SA_tol )];
    end
end

index = index(1:end);

%outputs
SA_out = SA(index);
FY_out = FY(index);
FZ_out = FZ(index);
P_out = P(index);
TSTC_out = TSTC(index);
IA_out = IA(index);
ET_out = ET(index);
NFY_out = NFY(index);
