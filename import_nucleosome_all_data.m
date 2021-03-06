function NucleosomeAll = import_nucleosome_all_data(filename, dataLines)
%IMPORTFILE Import data from a text file
%  NUCLEOSOMEALL = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  NUCLEOSOMEALL = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  NucleosomeAll = importfile("/Users/lawrimor/Downloads/200730_Nucleosome_All.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 30-Jul-2020 14:20:37

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 30);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["A", "D", "alpha", "area", "delta_area", "dist_err", "dist_to_53bp1", "dist_to_boundary", "exp_label", "frame", "mass", "mean", "particle", "peak", "phi", "r", "raw_data", "sig_raw", "sig_x", "sig_y", "sigx_to_sigraw", "sigy_to_sigraw", "sort_flag_53bp1", "sort_flag_boundary", "std", "traj_length", "x", "x_raw", "y", "y_raw"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "string", "double", "categorical", "double", "double", "double", "double", "double", "double", "double", "categorical", "double", "double", "double", "double", "double", "string", "string", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["dist_to_53bp1", "sort_flag_53bp1", "sort_flag_boundary"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["dist_to_53bp1", "exp_label", "raw_data", "sort_flag_53bp1", "sort_flag_boundary"], "EmptyFieldRule", "auto");

% Import the data
NucleosomeAll = readtable(filename, opts);

end