include("scripts.jl")
Random.seed!(777)

function scr_ml(ARGS)

    # PARSING ARGUMENTS
    # id_scr = 0; n_distr = 10; n_bio = 3; n_sample = 20; mxcb = 100
    # inp = "../data/data_f_cs.csv"; out = "../res"
    id_scr = parse(Int32, ARGS[1])
    # number of biotracers
    id_replicate = parse(Int64, ARGS[2])
    # number of indiv used for building distribution
    n_distr = parse(Int64, ARGS[3])
    # number of sample used to test
    n_sample = parse(Int64, ARGS[4])
    # number of replicates
    n_rep = parse(Int64, ARGS[5])
    # maximum column combintaions
    mxcb = parse(Int64, ARGS[6])
    # path to data (pca-transformed or not)
    pca = (parse(Int32, ARGS[7]) == 1)
    println(pca)
    inp = pca ? "../data/data_f_pca.csv" : "../data/data_f_cs.csv"
    if pca
        println("Only one set of columns will be used")
    end

    Random.seed!(randperm(10000)[id_replicate])

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

    out = "res_all123_"*ARGS[2]*"_"*ARGS[7]*"/"
    mkdir(out)
    println("Data loaded ✓")

    for l in 1:3
        # NB n_bio used as n_bio_max
        scr_nbio_v2(df, l, n_distr, n_sample, mxcb, n_rep, n_tot_bio,
        n_tot_indiv, n_reg, aug, sig, out, pca)
    end

    return 0
end


println("All scripts loaded ✓")
scr_ml(ARGS)

