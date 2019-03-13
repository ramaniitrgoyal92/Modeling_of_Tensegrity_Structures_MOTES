function [ output_structure ] = tenseg_field_init( field,input_structure,output_structure,default,print )
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% [ output_structure ] = TENSEG_FIELD_INIT( field,input_structure,output_structure,default,print )
% initialize the value for a field of a given structure by either using
% the value of the field variable in the workspace or loading the specified
% default value. Default values are specified in tenseg_defaults.
%
% Inputs:
%	field: name of field being initialized. If field name exists as
%		variable in workspace, that value will be used. Otherwise, the
%		corresponding field within default will be used.
%	input_structure: structure in which to check for field
%	output_structure: name of structure in which field is being initialized
%	default: structure containing all default field/variable values. These
%		are specified in tenseg_defaults.
%	print (optional): 1 or 0 specifying whether or not notification text is
%		printed to command window
%
% Outputs:
%	output_structure: the resulting modified structure


% Handle optional inputs
if nargin == 3
	print = 1;
end

% Check if specified field name exists as a variable in the input struct
try
	val = input_structure.(field);
catch
	val = [];
end

% If variable exists, store it as the field value
if ~isempty(val)
	output_structure.(field) = val;
	if print
		if numel(output_structure.(field))==1
			disp([field ': ' num2str(output_structure.(field))])
		else
			disp([field ':'])
			disp(output_structure.(field));
		end
	end
	
% If variable doesn't exist, store the default value
else
	output_structure.(field) = default.(field);
	if print
		if numel(output_structure.(field))==1
			disp([field ': ' num2str(output_structure.(field)) ' DEFAULT'])
		else
			disp([field ':' 'DEFAULT'])
			disp(output_structure.(field));
		end
	end
end

end

