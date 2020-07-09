include("scripts.jl")


function scr_ml(ARGS)
    # id_scr = 0; n_distr = 10; n_bio = 3; n_sample = 20; mxcb = 100
    # inp = "../data/data_f_cs.csv"; out = "../res"
    id_scr = parse(Int32, ARGS[1])
    # number of indiv used for building distribution
    n_distr = parse(Int64, ARGS[2])
    # number of biotracers
    n_bio = parse(Int64, ARGS[3])
    # number of sample used to test
    n_sample = parse(Int64, ARGS[4])
    # number of replicates
    n_rep = parse(Int64, ARGS[5])
    # maximum column combintaions
    mxcb = parse(Int64, ARGS[6])
    nrep = parse(Int64, ARGS[7])
    # path to data (pca-transformed or not)
    inp =  ARGS[8]
    # path to res folder
    out = ARGS[9]
    # number of regions considered
    n_reg = 3
    n_tot_indiv = 30
    n_tot_bio = 17
    # data augmentation parameters
    aug = Int64(1e3)
    sig = .1

    df = get_data(false, inp)

    println("Data loaded ✓")

    if id_scr == 1
        # NB n_distr used as n_distr_max
        scr_ndistr(df, n_bio, n_distr, n_sample, mxcb, n_rep, n_tot_bio, n_tot_indiv,
            n_reg, aug, sig, out)
    elseif id_scr == 2
        scr_nsample(df, n_bio, n_distr, n_sample, mxcb, n_rep, n_tot_bio,
        n_tot_indiv, n_reg, aug, sig, out)
    elseif id_scr == 3
        scr_noise(df, n_bio, n_distr, n_sample, mxcb, n_rep, n_tot_bio,
        n_tot_indiv, n_reg, aug, sig, out)
    else
        println("cool")
        # NB n_bio used as n_bio_max
        scr_nbio(df, n_bio, n_distr, n_sample, mxcb, n_rep, n_tot_bio,
        n_tot_indiv, n_reg, aug, sig, out)
    end

    return 0
end


println("All scripts loaded ✓")
scr_ml(ARGS)

# df = get_data(false, inp)
# inp = "../data/data_f_cs.csv"
# scr_train(df)


