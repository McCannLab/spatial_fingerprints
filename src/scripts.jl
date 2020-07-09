include("helpers.jl")


function scr_nbio(df, n_bio_max, n_distr, n_sample, mxcb, n_rep, n_tot_bio,
    n_tot_indiv, n_reg, aug, sig, path)
    for i in 1:n_bio_max
        println("------- n_bio = ", i, "------")
        cb_cols = sample_cols(i, n_tot_bio, mxcb)
        for j in 1:size(cb_cols, 1)
            println("  -- combin = ", j, "--")
            for k in 1:n_rep
                println("   -> replicate #", k)
                # all samples that will be used
                spl0 = samp_samples(n_distr + n_sample, n_tot_indiv, n_reg)
                # spl1 => sample in spl0 that'll be used in the training set
                spl1 = samp_samples(n_distr, n_distr + n_sample, n_reg)
                df_tmp = df[:, cb_cols[i]]
                # add origin (last column) back
                df_tmp.origin = df.origin
                fl = @sprintf "%sml_nbio_%02d_%03d_%03d.txt" path i j k
                # println(fl)
                res = simu_ml(df_tmp[spl0[spl1], :],
                    df_tmp[setdiff(spl0, spl0[spl1]), :], sig, fl, aug)
            end
        end
    end
end


function scr_ndistr(df, n_bio, n_distr_max, n_sample, mxcb, n_rep, n_tot_bio,
    n_tot_indiv, n_reg, aug, sig, path)
    for i in 1:n_distr_max
        println("------- n_distr = ", i, "------")
        cb_cols = sample_cols(n_bio, n_tot_bio, mxcb)
        for j in 1:size(cb_cols, 1)
            println("  -- combin = ", j, "--")
            for k in 1:n_rep
                println("   -> replicate #", k)
                # all samples that will be used
                spl0 = samp_samples(i + n_sample, n_tot_indiv, n_reg)
                # spl1 => sample in spl0 that'll be used in the training set
                spl1 = samp_samples(i, i + n_sample, n_reg)
                df_tmp = df[:, cb_cols[i]]
                df_tmp.origin = df.origin
                fl = @sprintf "%sml_ndistr_%02d_%03d_%03d.txt" path i j k
                res = simu_ml(df_tmp[spl0[spl1], :],
                    df_tmp[setdiff(spl0, spl0[spl1]), :], sig, fl, aug)
            end
        end
    end
end


function scr_nsample(df, n_bio, n_distr, n_sample_max, mxcb, n_rep, n_tot_bio,
    n_tot_indiv, n_reg, aug, sig, path)
    for i in 1:n_sample_max
        println("------- n_sample = ", i, "------")
        cb_cols = sample_cols(n_bio, n_tot_bio, mxcb)
        for j in 1:size(cb_cols, 1)
            println("  -- combin = ", j, "--")
            for k in 1:n_rep
                println("   -> replicate #", k)
                # all samples that will be used
                spl0 = samp_samples(n_distr + i, n_tot_indiv, n_reg)
                # spl1 => sample in spl0 that'll be used in the training set
                spl1 = samp_samples(n_distr, n_distr + i, n_reg)
                df_tmp = df[:, cb_cols[i]]
                df_tmp.origin = df.origin
                fl = @sprintf "%sml_nsamp_%02d_%03d_%03d.txt" path i j k
                res = simu_ml(df_tmp[spl0[spl1], :],
                    df_tmp[setdiff(spl0, spl0[spl1]), :], sig, fl, aug)
            end
        end
    end
end


function scr_noise(df, n_bio, n_distr, n_sample, mxcb, n_rep, n_tot_bio,
    n_tot_indiv, n_reg, aug, sig, path)
    noise = [exp10(k) for k in -4:.25:0]
    for (i, s) in enumerate(noise)
        println("- noise = ", i, "---------")
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
                df_tmp = df_tmp0[:, cb_cols[i]]
                df_tmp.origin = df.origin
                fl = @sprintf "%sml_noise_%02d_%03d_%03d.txt" path i j k
                res = simu_ml(df_tmp[spl0[spl1], :],
                    df_tmp[setdiff(spl0, spl0[spl1]), :], sig, fl, aug)
            end
        end
    end
end



function scr_train(df)
    sig = [exp10(k) for k in -4:.5:1]
    aug = 10*[5 10 50 100 500 1000 5000 10000 50000 100000]
    for (i, s) in enumerate(sig)
        for (j, a) in enumerate(aug)
        println("- train: i = ", i, "---- j = ", j)
        fl = @sprintf "../res_ml/train_%03d_%03d.txt" i j
        res = simu_ml(
            df[vcat(1:20, 31:50, 61:80), :],
            df[vcat(21:30, 51:60, 81:90), :],
            s, fl, a)
        end
    end
end
