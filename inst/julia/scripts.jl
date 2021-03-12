include("helpers.jl")

# Functions here have similar structure and thus similar arguments:
# details arguments here:
# - df : data frame include all bio-tracers values and tru origins

# increase the number of bio-tracers
function scr_nbio_v2(df, n_bio, n_distr, n_sample, mxcb, n_rep, n_tot_bio,
    n_tot_indiv, n_reg, aug, sig, path, pca = false)

    if pca
        cb_cols = [1:n_bio]
    else
        cb_cols = sample_cols(n_bio, n_tot_bio, mxcb)
    end

    # loop over bio-tracers combinations
    for j in 1:size(cb_cols, 1)
        println("  -- combin = ", j, "--")
        # loop over repetition
        for k in 1:n_rep
            # println("   -> replicate #", k)
            # all samples that will be used
            spl0 = samp_samples(n_distr + n_sample, n_tot_indiv, n_reg)
            # spl1 => sample in spl0 that'll be used in the training set
            spl1 = samp_samples(n_distr, n_distr + n_sample, n_reg)
            df_tmp = df[:, cb_cols[j]]
            # add origin (last column) back
            df_tmp.origin = df.origin
            fl = @sprintf "%sml_nbio_%02d_%03d_%05d.txt" path n_bio j k
            train = df_tmp[spl0[spl1], :]
            valid = df_tmp[setdiff(spl0, spl0[spl1]), :]
            res = simu_ml(train, valid, sig, fl, aug)
        end
    end
end

# increasing the number of bio-tracers
function scr_nbio(df, n_bio_max, n_distr, n_sample, mxcb, n_rep, n_tot_bio,
    n_tot_indiv, n_reg, aug, sig, path, pca = false)
    for i in 1:n_bio_max
        println("------- n_bio = ", i, "------")
        if pca
            cb_cols = [1:n_bio_max]
        else
            cb_cols = sample_cols(i, n_tot_bio, mxcb)
        end
        for j in 1:size(cb_cols, 1)
            println("  -- combin = ", j, "--")
            for k in 1:n_rep

                # println("   -> replicate #", k)
                # all samples that will be used
                spl0 = samp_samples(n_distr + n_sample, n_tot_indiv, n_reg)
                # spl1 => sample in spl0 that'll be used in the training set
                spl1 = samp_samples(n_distr, n_distr + n_sample, n_reg)
                df_tmp = df[:, cb_cols[j]]
                # add origin (last column) back
                df_tmp.origin = df.origin
                fl = @sprintf "%sml_nbio_%02d_%03d_%03d.txt" path i j k
                train = df_tmp[spl0[spl1], :]
                valid = df_tmp[setdiff(spl0, spl0[spl1]), :]
                res = simu_ml(train, valid, sig, fl, aug)
            end
        end
    end
end


# increasing the size of the training set
function scr_ndistr(df, n_bio, n_sample, mxcb, n_rep, n_tot_bio, n_tot_indiv,
    n_reg, aug, sig, id_rep)
    # set seed
    Random.seed!(randperm(1000)[id_rep])
    # loop over training set size
    for i in 5:25
        println("------- n_distr = ", i, "------")
        cb_cols = sample_cols(n_bio, n_tot_bio, mxcb)
        for j in 1:size(cb_cols, 1)
            println("  -- combin = ", j, "--")
            # loop over repetition
            for k in 1:n_rep
                println("   -> replicate #", k)
                # all samples that will be used
                spl0 = samp_samples(i + n_sample, n_tot_indiv, n_reg)
                # spl1 => sample in spl0 that'll be used in the training set
                spl1 = samp_samples(i, i + n_sample, n_reg)
                df_tmp = df[:, cb_cols[j]]
                df_tmp.origin = df.origin
                train = df_tmp[spl0[spl1], :]
                valid = df_tmp[setdiff(spl0, spl0[spl1]), :]
                fl = @sprintf "ml_ndistr_%02d_%02d.txt" n_bio id_rep
                res = simu_ml2(train, valid, sig, fl, i, j, k, aug)
            end
        end
    end
end


# increasing the number of sample (not used anymore as we used nbio instead)
function scr_nsample(df, n_bio, n_distr, n_sample_max, mxcb, n_rep, n_tot_bio,
    n_tot_indiv, n_reg, aug, sig, path, pca = false)
    for i in 1:n_sample_max
        println("------- n_sample = ", i, "------")
        if pca
            cb_cols = [1:n_bio_max]
        else
            cb_cols = sample_cols(n_bio, n_tot_bio, mxcb)
        end
        for j in 1:size(cb_cols, 1)
            println("  -- combin = ", j, "--")
            for k in 1:n_rep
                println("   -> replicate #", k)
                # all samples that will be used
                spl0 = samp_samples(n_distr + i, n_tot_indiv, n_reg)
                # spl1 => sample in spl0 that'll be used in the training set
                spl1 = samp_samples(n_distr, n_distr + i, n_reg)
                df_tmp = df[:, cb_cols[j]]
                df_tmp.origin = df.origin
                fl = @sprintf "%sml_nsamp_%02d_%03d_%03d.txt" path i j k
                res = simu_ml(df_tmp[spl0[spl1], :],
                    df_tmp[setdiff(spl0, spl0[spl1]), :], sig, fl, aug)
            end
        end
    end
end


# increase noise level
function scr_noise(df, n_bio, n_distr, n_sample, mxcb, n_rep, n_tot_bio,
    n_tot_indiv, n_reg, aug, sig, id_rep)
    # set seed
    Random.seed!(randperm(1000)[id_rep])
    # noise sequence
    noise = [exp10(k) for k in -4:.2:1]
    for (i, s) in enumerate(noise)
        println("- noise = ", s, "---------")
        cb_cols = sample_cols(n_bio, n_tot_bio, mxcb)
        for j in 1:size(cb_cols, 1)
            println("  -- combin = ", j, "--")
            for k in 1:n_rep
                println("   -> replicate #", k)
                df_tmp0 = add_noise_df(df, s)
                # all samples that will be used
                spl0 = samp_samples(n_distr + n_sample, n_tot_indiv, n_reg)
                # spl1 => sample in spl0 that'll be used in the training set
                spl1 = samp_samples(n_distr, n_distr + n_sample, n_reg)
                df_tmp = df_tmp0[:, cb_cols[j]]
                df_tmp.origin = df.origin
                fl = @sprintf "ml_noise_%02d_%02d.txt" n_bio id_rep
                res = simu_ml2(df_tmp[spl0[spl1], :],
                    df_tmp[setdiff(spl0, spl0[spl1]), :], sig, fl, s, j, k, aug)
            end
        end
    end
end


# testing the influence if noise addition and data augmentation
function scr_train(df, col, fl)
    sig = [0 0.00001 0.0001 0.001 0.005 0.01 0.05 0.1 .2 0.5 1 2]
    aug = [1 2 5 10 20 50 100 200 500 1000 2000 5000]
    for (i, s) in enumerate(sig)
        for (j, a) in enumerate(aug)
        println("- train: i = ", i, "---- j = ", j)
        res = simu_ml2(
            df[vcat(1:20, 31:50, 61:80), col],
            df[vcat(21:30, 51:60, 81:90), col],
            s, fl, size(col, 1), s, a, a)
        end
    end
end

# df = get_data(false, "../data/data_f_cs.csv")
# scr_train(df, 1:18, "res_train_all.txt")
# scr_train(df, vcat(1:3, 18), "res_train_si.txt")
# scr_train(df, vcat(4:6, 18), "res_train_fa.txt")
