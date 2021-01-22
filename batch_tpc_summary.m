function S_tpc_summary = batch_tpc_summary(S_flat, pixel_size_um, min_overlap)

fnames = fieldnames(S_flat);

for f = 1:numel(fnames)
    S_tpc_summary.(fnames{f}) = bin_tpc(S_flat.(fnames{f}), pixel_size_um, min_overlap);
end
S_tpc_summary.min_overlap = min_overlap;
S_tpc_summary.pixel_size_um = pixel_size_um; 