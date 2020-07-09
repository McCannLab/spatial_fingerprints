using Parameters, CSV, DataFrames, DataFramesMeta, Query
using Statistics, Random, Distributions, Combinatorics, Flux, Tracker
using Printf, JLD, Plots




# add Parameters CSV DataFrames Statistics Random Distributions
# add DataFramesMeta Flux Plots

myscale(A) = (A .- mean(A))./std(A)

function get_data(pca = false, path = "data/data_f_pca.csv")
    if (pca)
        df = CSV.File(path) |> DataFrame!
    else
        df = CSV.File(path) |> DataFrame!
    end
    select!(df, Not(1:2))
    df.origin = repeat(1:3, inner = 30)
    return df
end


function simu_accuracy(df_test, mod, nreg::Int64 = 3)
    df_tmp = select(df_test, Not(:origin))
    out = zeros(Float32, nreg, nreg)
    for i in 1:size(df_tmp, 1)
        new = [df_tmp[i, k] for k in 1:size(df_tmp, 2)]
        res = mod(new)
        # Tracker.data ==> otherwise we have tracker array that cannot be
        # use with "base array"
        out[df_test.origin[i], :] += Tracker.data(res)
    end
    return out
end
# df4 = df2[:,1:3]
# df4.origin = df2.origin
# simu_accuracy(df4, m)

# returns a vector of n indiv sample among ndiv for nreg regions
function samp_samples(n::Int64, nindiv::Int64 = 30, nreg::Int64 = 3)
    out = []
    for i in 1:nreg
        append!(out, randperm(nindiv)[1:n] .+ (i - 1) * nindiv)
    end
    return out
end

function write_res(df_test, mod, file)
    df_tmp = select(df_test, Not(:origin))
    io = open(file, "w")
    for i in 1:size(df_tmp, 1)
        new = [df_tmp[i, k] for k in 1:size(df_tmp, 2)]
        res = Tracker.data(mod(new))
        for j in 1:size(res, 1)
            write_res_unit(i, j, df_test.origin[i], res[j], io)
        end
    end
    close(io)
    return nothing
end

function write_res_unit(id_sampl, reg_hyp, reg_true, val, io)
    println(io, id_sampl, " ", reg_hyp, " ", reg_true, " ", val)
    return nothing
end

function add_noise(val, σ)
    val .+ σ * randn(size(val, 1))
end

function add_noise_df(df, σ)
    DataFrame([add_noise(col, σ) for col = eachcol(df)])
end



# taining set with data augmentation (white noise)
function simu_train(df_train, sig, nrep::Int64 = 10, nreg::Int64 = 3)
    df_tmp = select(df_train, Not(:origin))
    train_data = []
    for r in 1:nrep
        for j in 1:size(df_train, 1)
            x = [df_tmp[j, k] for k = 1:size(df_tmp, 2)]
            x += sig * randn(size(df_tmp, 2))
            y = fill(0.0, nreg)
            y[df_train.origin[j]] = 1.0
            push!(train_data, (x, y))
        end
    end
    return train_data
end
# simu_train(df2, .1, 2)
# df3 = select(df2, Not(:origin))
# simu_accuracym


# columns sampler
function sample_cols(nbio, mxbio, mxcb)
    tmp = collect(combinations(Array(1:mxbio), nbio))
    if (size(tmp, 1) <= mxcb)
        out = tmp
    else
        out = tmp[randperm(size(tmp, 1))[1:mxcb]]
    end
    out
end


function simu_ml(df_train, df_test, sig, file::String, aug = 1e5,
    nreg::Int64 = 3)
    ##
    ndim = size(df_train, 2) - 1
    nhid = round(Int, (ndim + nreg) / 2)
    ##
    m = Chain(Dense(ndim, nhid, σ), Dense(nhid, nreg), softmax)
    # loss(x, y) = Flux.mse(m(x), y)
    loss(x, y) = Flux.crossentropy(m(x), y)
    #
    ##
    print("Generating data set ")
    train_data = simu_train(df_train, sig, aug)
    println("✓")
    ##
    print("Training ")
    Flux.train!(loss, Flux.params(m), train_data, ADAM())
    println("✓")
    ##
    write_res(df_test, m, file)
    ##
    return nothing
end





