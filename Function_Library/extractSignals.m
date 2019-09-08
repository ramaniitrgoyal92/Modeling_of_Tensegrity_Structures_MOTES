function signals = extractSignals(odefun, names)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
%
%EXTRACTSIGNALS Get log signals from ODE output.
%   SIGNALS = EXTRACTSIGNALS(@ODEFUN, NAMES) evaluates the ODEFUN to get the
%       auxiliary output signal logs, and assigns them to a structure with
%       the signal names listed in NAMES.
%
%   Note that it is assumed that the NAMES are provided in the same order
%   as the signals in the log output. If more NAMES are provided than
%   signals in the log cell produced by ODEFUN, the additional names will
%   be ignored.
%

% Get ODE output. Assumes that there are three ODEFUN inputs: time,
% state, and parameters.
    [~, logsout] = odefun([], [], []);

    % Get number of signals and number of named signals.
    nSignals = length(logsout);
    nNames = length(names);

    % If we find no log signals, throw an error.
    if ~nSignals
        error('extractSignals:NoSignals', ...
              'No log signals found. Re-run ODE solver!')
    end

    % Create output structure to hold log variables.
    signals = struct();

    for idx = 1 : nNames
        if idx <= nSignals
            signals.(names{idx}) = logsout{idx};
        end
    end
end
