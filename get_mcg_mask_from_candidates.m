function mask = get_mcg_mask_from_candidates(candidates_scg, candidate_id)
    mask = ismember(candidates_scg.superpixels, candidates_scg.labels{candidate_id});

