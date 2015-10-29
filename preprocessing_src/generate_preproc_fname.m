function fname = generate_preproc_fname(dataset, video_name, prepr_params)


prepr_params.dataset = dataset;
prepr_params.video = video_name;

fname = regexprep(sprintf('ppvid_%s', buildStringFromStruct(prepr_params, '__')), '\.', '_' );
