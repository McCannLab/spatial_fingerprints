include("scripts.jl")


function scr_ml(ARGS)

    # PARSING ARGUMENTS
    id_scr = parse(Int32, ARGS[1])
    # number of biotracers
    n_bio = parse(Int64, ARGS[2])
    # number of indiv used for building distribution
    n_distr = parse(Int64, ARGS[3]) # NB: not used here
    # number of sample used to test
    n_sample = parse(Int64, ARGS[4])
    # number of replicates
    n_rep = parse(Int64, ARGS[5])
    # maximum column combintaions
    mxcb = parse(Int64, ARGS[6])
    # path to data (pca-transformed or not)
    id_rep = parse(Int64, ARGS[7])
    # path to data
    inp="../data/data_f_cs.csv"
    # FIXED PARAMETERS
    # path to res folder
    # number of regions considered
    n_reg = 3
    n_tot_indiv = 30
    n_tot_bio = 17
    # data augmentation parameters
    aug = Int64(1e3)
    sig = .01

    df = get_data(false, inp)

    out = "./"

    println("Data loaded ✓")

    if id_scr == 1
        # NB n_distr is ignored
        scr_ndistr(df, n_bio, n_sample, mxcb, n_rep, n_tot_bio, n_tot_indiv,
            n_reg, aug, sig, id_rep)
    elseif id_scr == 2
        scr_nsample(df, n_bio, n_distr, n_sample, mxcb, n_rep, n_tot_bio,
        n_tot_indiv, n_reg, aug, sig, out, pca)
    elseif id_scr == 3
        scr_noise(df, n_bio, n_distr, n_sample, mxcb, n_rep, n_tot_bio,
        n_tot_indiv, n_reg, aug, sig, id_rep)
    else
        # NB n_bio used as n_bio_max
        scr_nbio_v2(df, n_bio, n_distr, n_sample, mxcb, n_rep, n_tot_bio,
        n_tot_indiv, n_reg, aug, sig, out, pca)
    end

    return 0
end


println("All scripts loaded ✓")
scr_ml(ARGS)
